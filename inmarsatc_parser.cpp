#include <inmarsatc_parser.h>

namespace inmarsatc {
    namespace frameParser {
        //START CLASS PacketDecoder
        PacketDecoder::packetDecoder_result PacketDecoder::basicDecode(decoder::Decoder::decoder_result inputFrame, int* pos) {
            packetDecoder_result ret;
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_NONE;
            // As a safe precaution, in case we cannot correctly ascertain the packet length,
            // setting the length this way will allow us to discard the whole remaining frame
            // and avoid an infinite loop.
            ret.packetLength = 640 - *pos;
            // Packet descriptor
            ret.packetDescriptor = inputFrame.decodedFrame[*pos];
            /// There are 2 variations of packet descriptor
            /// Short packet descriptor and Medium packet descriptor
            /// They give us the packet lengths
            if (ret.packetDescriptor >> 7 == 0) {
                /// Short packet descriptor
                /// The packet length including CRC does not include byte 0, we add 1
                ret.packetLength = (ret.packetDescriptor & 0x0F) + 1;
            } else if (ret.packetDescriptor >> 6 == 0x02) {
                /// Medium packet descriptor
                /// The packet length including CRC does not include the first 2 bytes, we add 2
                ret.packetLength = inputFrame.decodedFrame[*pos + 1] + 2;
            }
            // At this stage we do not know for sure if the CRC is correct.
            /// We compute the 2-byte CRC and compare with the packet 2-byte CRC
            /// The packet 2-byte CRC position is given by the packet descriptor
            int packetCRC = (inputFrame.decodedFrame[*pos + ret.packetLength - 2] << 8) | inputFrame.decodedFrame[*pos + ret.packetLength - 1];
            int computedCRC = computeCRC(inputFrame.decodedFrame, *pos, ret.packetLength);
            //Added a check to zero.
            //The BD-BE packet content is a packet that does not have a CRC.
            //This is a workaround to avoid computing the CRC for the content.
            //This needs be fixed by actually computing the CRC. TODO
            ret.isCrc = packetCRC == 0 || packetCRC == computedCRC;
            ret.frameNumber = inputFrame.frameNumber;
            ret.packetVars.insert(std::pair<std::string, std::string>("packetDescriptorText", getDescriptorAsText(ret.packetDescriptor)));
            ret.timestamp = inputFrame.timestamp;
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_27(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_COMPLETE;
            int mesId = inputFrame.decodedFrame[*pos + 1] << 16 | inputFrame.decodedFrame[*pos + 1 + 1] << 8 | inputFrame.decodedFrame[*pos + 1 + 2];
            int sat = inputFrame.decodedFrame[*pos + 4] >> 6 & 0x03;
            std::string satName = getSatName(sat);
            int lesId = inputFrame.decodedFrame[*pos + 4] & 0x3F;
            std::string lesName = getLesName(sat, lesId);
            int logicalChannelNo = inputFrame.decodedFrame[*pos + 5];
            ret.packetVars.insert(std::pair<std::string, std::string>("mesId", std::to_string(mesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("sat", std::to_string(sat)));
            ret.packetVars.insert(std::pair<std::string, std::string>("satName", satName));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesId", std::to_string(lesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesName", lesName));
            ret.packetVars.insert(std::pair<std::string, std::string>("logicalChannelNo", std::to_string(logicalChannelNo)));
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_2A(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_PARTIAL;
            int mesId = inputFrame.decodedFrame[*pos + 1] << 16 | inputFrame.decodedFrame[*pos + 1 + 1] << 8 | inputFrame.decodedFrame[*pos + 1 + 2];
            int sat = inputFrame.decodedFrame[*pos + 4] >> 6 & 0x03;
            std::string satName = getSatName(sat);
            int lesId = inputFrame.decodedFrame[*pos + 4] & 0x3F;
            std::string lesName = getLesName(sat, lesId);
            int logicalChannelNo = inputFrame.decodedFrame[*pos + 5];
            std::ostringstream os;
            for(int i = 0; i < 3; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << (uint16_t)inputFrame.decodedFrame[*pos + 6 + i];
            }
            std::string unknown1Hex = os.str();
            ret.packetVars.insert(std::pair<std::string, std::string>("mesId", std::to_string(mesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("sat", std::to_string(sat)));
            ret.packetVars.insert(std::pair<std::string, std::string>("satName", satName));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesId", std::to_string(lesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesName", lesName));
            ret.packetVars.insert(std::pair<std::string, std::string>("logicalChannelNo", std::to_string(logicalChannelNo)));
            ret.packetVars.insert(std::pair<std::string, std::string>("unknown1Hex", unknown1Hex));
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_08(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_COMPLETE;
            int sat = inputFrame.decodedFrame[*pos + 1] >> 6 & 0x03;
            std::string satName = getSatName(sat);
            int lesId = inputFrame.decodedFrame[*pos + 1] & 0x3F;
            std::string lesName = getLesName(sat, lesId);
            int logicalChannelNo = inputFrame.decodedFrame[*pos + 2];
            double uplinkChannelMhz = ((inputFrame.decodedFrame[*pos + 3] << 8 | (uint16_t)inputFrame.decodedFrame[*pos + 3 + 1]) - 6000) * 0.0025 + 1626.5;
            ret.packetVars.insert(std::pair<std::string, std::string>("sat", std::to_string(sat)));
            ret.packetVars.insert(std::pair<std::string, std::string>("satName", satName));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesId", std::to_string(lesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesName", lesName));
            ret.packetVars.insert(std::pair<std::string, std::string>("logicalChannelNo", std::to_string(logicalChannelNo)));
            ret.packetVars.insert(std::pair<std::string, std::string>("uplinkChannelMhz", std::to_string(uplinkChannelMhz)));
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_6C(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_COMPLETE;
            //The second byte in 6C packets is the services byte.
            //Identical to the first byte of the services bytes of the 7D packet.
            int is8 = inputFrame.decodedFrame[*pos + 1];
            std::string services = getServices_short(is8);
            double uplinkChannelMhz = ((inputFrame.decodedFrame[*pos + 2] << 8 | (uint16_t)inputFrame.decodedFrame[*pos + 2 + 1]) - 6000) * 0.0025 + 1626.5;
            int tdmslots_int[28];
            int j = 0;
            for (int i = 0; i < 28; i += 4) {
                tdmslots_int[i] = inputFrame.decodedFrame[*pos + 4 + j] >> 6;
                tdmslots_int[i + 1] = inputFrame.decodedFrame[*pos + 4 + j] >> 4 & 3;
                tdmslots_int[i + 2] = inputFrame.decodedFrame[*pos + 4 + j] >> 2 & 3;
                tdmslots_int[i + 3] = inputFrame.decodedFrame[*pos + 4 + j] & 3;
                j++;
            }
            std::string tdmSlots;
            for(int i = 0; i < 28; i++) {
                tdmSlots += std::to_string(i) + ": " + std::to_string(tdmslots_int[i]) + "\n";
            }
            ret.packetVars.insert(std::pair<std::string, std::string>("services", services));
            ret.packetVars.insert(std::pair<std::string, std::string>("uplinkChannelMhz", std::to_string(uplinkChannelMhz)));
            ret.packetVars.insert(std::pair<std::string, std::string>("tdmSlots", tdmSlots));
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_7D(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_COMPLETE;
            int networkVersion = inputFrame.decodedFrame[*pos + 1];
            ret.frameNumber = inputFrame.decodedFrame[*pos + 2] << 8 | inputFrame.decodedFrame[*pos + 3];
            double timestamp_seconds_d = ret.frameNumber * 8.64;
            int timestamp_hours = floor(timestamp_seconds_d/3600.0);
            int timestamp_min = floor((((int)timestamp_seconds_d)%3600)/60.0);
            int timestamp_sec = ((int)timestamp_seconds_d)%60;
            int timestamp_msec = (((int)(timestamp_seconds_d*1000))%1000);
            std::string timestamp_str = std::to_string(timestamp_hours) + ":" + std::to_string(timestamp_min) + ":" + std::to_string(timestamp_sec) + "." + std::to_string(timestamp_msec);
            int signallingChannel = inputFrame.decodedFrame[*pos + 4] >> 2;
            int count = (inputFrame.decodedFrame[*pos + 5] >> 4 & 0x0F) * 0x02;
            int channelType = inputFrame.decodedFrame[*pos + 6] >> 0x05;
            std::string channelTypeName;
            switch (channelType) {
                case 1:
                    channelTypeName = "NCS";
                case 2:
                    channelTypeName = "LES TDM";
                case 3:
                    channelTypeName = "Joint NCS and TDM";
                case 4:
                    channelTypeName = "ST-BY NCS";
                default:
                    channelTypeName = "Reserved";
            }
            int local = inputFrame.decodedFrame[*pos + 6] >> 2 & 0x07;
            int sat = inputFrame.decodedFrame[*pos + 7] >> 6 & 0x03;
            std::string satName = getSatName(sat);
            int lesId = inputFrame.decodedFrame[*pos + 7] & 0x3F;
            std::string lesName = getLesName(sat, lesId);
            uint8_t status_b = inputFrame.decodedFrame[*pos + 8];
            std::string status;
            status += "Bauds600: " + std::to_string((status_b & 0x80) >> 7 == 1) + "\n";
            status += "Operational: " + std::to_string((status_b & 0x40) >> 6 == 1) + "\n";
            status += "InService: " + std::to_string((status_b & 0x20) >> 5 == 1) + "\n";
            status += "Clear: " + std::to_string((status_b & 0x10) >> 4 == 1) + "\n";
            status += "LinksOpen: " + std::to_string((status_b & 0x08) >> 3 == 1);
            int iss = inputFrame.decodedFrame[*pos + 9] << 8 | inputFrame.decodedFrame[*pos + 10];
            std::string services = getServices(iss);
            int randomInterval = inputFrame.decodedFrame[*pos + 11];
            ret.packetVars.insert(std::pair<std::string, std::string>("networkVersion", std::to_string(networkVersion)));
            ret.packetVars.insert(std::pair<std::string, std::string>("sat", std::to_string(sat)));
            ret.packetVars.insert(std::pair<std::string, std::string>("satName", satName));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesId", std::to_string(lesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesName", lesName));
            ret.packetVars.insert(std::pair<std::string, std::string>("signallingChannel", std::to_string(signallingChannel)));
            ret.packetVars.insert(std::pair<std::string, std::string>("count", std::to_string(count)));
            ret.packetVars.insert(std::pair<std::string, std::string>("channelType", std::to_string(channelType)));
            ret.packetVars.insert(std::pair<std::string, std::string>("channelTypeName", channelTypeName));
            ret.packetVars.insert(std::pair<std::string, std::string>("local", std::to_string(local)));
            ret.packetVars.insert(std::pair<std::string, std::string>("status", status));
            ret.packetVars.insert(std::pair<std::string, std::string>("services", services));
            ret.packetVars.insert(std::pair<std::string, std::string>("randomInterval", std::to_string(randomInterval)));
            ret.packetVars.insert(std::pair<std::string, std::string>("timestamp_str", timestamp_str));
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_81(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_PARTIAL;
            int mesId = inputFrame.decodedFrame[*pos + 2] << 16 | inputFrame.decodedFrame[*pos + 2 + 1] << 8 | inputFrame.decodedFrame[*pos + 1 + 2];
            int sat = inputFrame.decodedFrame[*pos + 5] >> 6 & 0x03;
            std::string satName = getSatName(sat);
            int lesId = inputFrame.decodedFrame[*pos + 5] & 0x3F;
            std::string lesName = getLesName(sat, lesId);
            int logicalChannelNo = inputFrame.decodedFrame[*pos + 9];
            double downlinkChannelMhz = ((inputFrame.decodedFrame[*pos + 6] << 8 | inputFrame.decodedFrame[*pos + 6 + 1]) - 8000) * 0.0025 + 15305;
            int presentation = inputFrame.decodedFrame[*pos + 14];
            std::ostringstream os;
            for(int i = 0; i < 1; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << (uint16_t)inputFrame.decodedFrame[*pos + 8 + i];
            }
            std::string unknown1Hex = os.str();
            os.clear();
            for(int i = 0; i < 4; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << (uint16_t)inputFrame.decodedFrame[*pos + 10 + i];
            }
            std::string unknown2Hex = os.str();
            os.clear();
            for(int i = 0; i < 2; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << (uint16_t)inputFrame.decodedFrame[*pos + 15 + i];
            }
            std::string unknown3Hex = os.str();
            ret.packetVars.insert(std::pair<std::string, std::string>("mesId", std::to_string(mesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("sat", std::to_string(sat)));
            ret.packetVars.insert(std::pair<std::string, std::string>("satName", satName));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesId", std::to_string(lesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesName", lesName));
            ret.packetVars.insert(std::pair<std::string, std::string>("logicalChannelNo", std::to_string(logicalChannelNo)));
            ret.packetVars.insert(std::pair<std::string, std::string>("downlinkChannelMhz", std::to_string(downlinkChannelMhz)));
            ret.packetVars.insert(std::pair<std::string, std::string>("presentation", std::to_string(presentation)));
            ret.packetVars.insert(std::pair<std::string, std::string>("unknown1Hex", unknown1Hex));
            ret.packetVars.insert(std::pair<std::string, std::string>("unknown2Hex", unknown2Hex));
            ret.packetVars.insert(std::pair<std::string, std::string>("unknown3Hex", unknown3Hex));
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_83(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_PARTIAL;
            int mesId = inputFrame.decodedFrame[*pos + 2] << 16 | inputFrame.decodedFrame[*pos + 2 + 1] << 8 | inputFrame.decodedFrame[*pos + 2 + 2];
            int sat = inputFrame.decodedFrame[*pos + 5] >> 6 & 0x03;
            std::string satName = getSatName(sat);
            int lesId = inputFrame.decodedFrame[*pos + 5] & 0x3F;
            std::string lesName = getLesName(sat, lesId);
            uint8_t status_bits = inputFrame.decodedFrame[*pos + 6];
            int logicalChannelNo = inputFrame.decodedFrame[*pos + 7];
            int frameLength = inputFrame.decodedFrame[*pos + 8];
            int duration = inputFrame.decodedFrame[*pos + 9];
            double downlinkChannelMhz = ((inputFrame.decodedFrame[*pos + 10] << 8 | inputFrame.decodedFrame[*pos + 10 + 1]) - 8000) * 0.0025 + 1530.5;
            double uplinkChannelMhz = ((inputFrame.decodedFrame[*pos + 12] << 8 | inputFrame.decodedFrame[*pos + 12 + 1]) - 6000) * 0.0025 + 1626.5;
            int frameOffset = inputFrame.decodedFrame[*pos + 14];
            uint8_t packetDescriptor1 = inputFrame.decodedFrame[*pos + 15];
            ret.packetVars.insert(std::pair<std::string, std::string>("mesId", std::to_string(mesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("sat", std::to_string(sat)));
            ret.packetVars.insert(std::pair<std::string, std::string>("satName", satName));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesId", std::to_string(lesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesName", lesName));
            ret.packetVars.insert(std::pair<std::string, std::string>("status_bits", std::to_string(status_bits)));
            ret.packetVars.insert(std::pair<std::string, std::string>("logicalChannelNo", std::to_string(logicalChannelNo)));
            ret.packetVars.insert(std::pair<std::string, std::string>("frameLength", std::to_string(frameLength)));
            ret.packetVars.insert(std::pair<std::string, std::string>("duration", std::to_string(duration)));
            ret.packetVars.insert(std::pair<std::string, std::string>("downlinkChannelMhz", std::to_string(downlinkChannelMhz)));
            ret.packetVars.insert(std::pair<std::string, std::string>("uplinkChannelMhz", std::to_string(uplinkChannelMhz)));
            ret.packetVars.insert(std::pair<std::string, std::string>("frameOffset", std::to_string(frameOffset)));
            ret.packetVars.insert(std::pair<std::string, std::string>("packetDescriptor1", std::to_string(packetDescriptor1)));
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_91(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            return ret; //not implemented yet
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_92(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_COMPLETE;
            int loginAckLength = inputFrame.decodedFrame[*pos + 1];
            std::ostringstream os;
            for(int i = 0; i < 3; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << (uint16_t)inputFrame.decodedFrame[*pos + 2 + i];
            }
            std::string les = os.str();
            double downlinkChannelMhz = ((inputFrame.decodedFrame[*pos + 5] << 8 | inputFrame.decodedFrame[*pos + 5 + 1]) - 8000) * 0.0025 + 1530.5;
            os.clear();
            for(int i = 0; i < 1; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << (uint16_t)inputFrame.decodedFrame[*pos + 7 + i];
            }
            std::string stationStartHex = os.str();
            ret.packetVars.insert(std::pair<std::string, std::string>("loginAckLength", std::to_string(loginAckLength)));
            ret.packetVars.insert(std::pair<std::string, std::string>("downlinkChannelMhz", std::to_string(downlinkChannelMhz)));
            ret.packetVars.insert(std::pair<std::string, std::string>("les", les));
            ret.packetVars.insert(std::pair<std::string, std::string>("stationStartHex", stationStartHex));
            if(loginAckLength > 7) {
                //stations
                int stationCount = inputFrame.decodedFrame[8];
                std::string stations = getStations(inputFrame.decodedFrame, stationCount, *pos + 9);
                ret.packetVars.insert(std::pair<std::string, std::string>("stationCount", std::to_string(stationCount)));
                ret.packetVars.insert(std::pair<std::string, std::string>("stations", stations));
            }
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_9A(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            return ret; //not implemented yet
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_A0(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            return ret; //not implemented yet
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_A3(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_PARTIAL;
            int mesId = inputFrame.decodedFrame[*pos + 2] << 16 | inputFrame.decodedFrame[*pos + 2 + 1] << 8 | inputFrame.decodedFrame[*pos + 2 + 2];
            int sat = inputFrame.decodedFrame[*pos + 5] >> 6 & 0x03;
            std::string satName = getSatName(sat);
            int lesId = inputFrame.decodedFrame[*pos + 5] & 0x3F;
            std::string lesName = getLesName(sat, lesId);
            std::ostringstream os;
            if(ret.packetLength >= 38) {
                int j = *pos + 13;
                std::string shortMessage;
                for(int i = 0; j < *pos + ret.packetLength - 2; i++) {
                    shortMessage += (char)inputFrame.decodedFrame[j] & 0x7F; //x-IA5 encoding
                    j++;
                }
                ret.packetVars.insert(std::pair<std::string, std::string>("shortMessage", shortMessage));
                for(int i = 0; i < 6; i++) {
                    os << std::setfill('0') << std::setw(2) << std::right << std::hex << inputFrame.decodedFrame[*pos + 6 + i];
                }
            } else {
                for(int i = 0; i < ret.packetLength - 6; i++) {
                    os << std::setfill('0') << std::setw(2) << std::right << std::hex << (uint16_t)inputFrame.decodedFrame[*pos + 6 + i];
                }
            }
            std::string unknown1Hex = os.str();
            ret.packetVars.insert(std::pair<std::string, std::string>("mesId", std::to_string(mesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("sat", std::to_string(sat)));
            ret.packetVars.insert(std::pair<std::string, std::string>("satName", satName));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesId", std::to_string(lesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesName", lesName));
            ret.packetVars.insert(std::pair<std::string, std::string>("unknown1Hex", unknown1Hex));
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_A8(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_PARTIAL;
            int mesId = inputFrame.decodedFrame[*pos + 2] << 16 | inputFrame.decodedFrame[*pos + 2 + 1] << 8 | inputFrame.decodedFrame[*pos + 2 + 2];
            int sat = inputFrame.decodedFrame[*pos + 5] >> 6 & 0x03;
            std::string satName = getSatName(sat);
            int lesId = inputFrame.decodedFrame[*pos + 5] & 0x3F;
            std::string lesName = getLesName(sat, lesId);
            std::ostringstream os;
            for(int i = 0; i < 3; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << (uint16_t)inputFrame.decodedFrame[*pos + 6 + i];
            }
            std::string unknown1Hex = os.str();
            int shortMessageLength = inputFrame.decodedFrame[*pos + 9];
            os.clear();
            for(int i = 0; i < 1; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << (uint16_t)inputFrame.decodedFrame[*pos + 10 + i];
            }
            std::string unknown2Hex = os.str();
            if(shortMessageLength > 2) {
                std::string shortMessage;
                int j = *pos + 11;
                for(int i = 0; j < *pos + ret.packetLength - 2; i++) {
                    shortMessage += (char)inputFrame.decodedFrame[j] & 0x7F; //x-IA5 encoding
                    j++;
                }
                ret.packetVars.insert(std::pair<std::string, std::string>("shortMessage", shortMessage));
            }
            ret.packetVars.insert(std::pair<std::string, std::string>("mesId", std::to_string(mesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("sat", std::to_string(sat)));
            ret.packetVars.insert(std::pair<std::string, std::string>("satName", satName));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesId", std::to_string(lesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesName", lesName));
            ret.packetVars.insert(std::pair<std::string, std::string>("unknown1Hex", unknown1Hex));
            ret.packetVars.insert(std::pair<std::string, std::string>("unknown2Hex", unknown2Hex));
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_AA(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_COMPLETE;
            int sat = inputFrame.decodedFrame[*pos + 2] >> 6 & 0x03;
            std::string satName = getSatName(sat);
            int lesId = inputFrame.decodedFrame[*pos + 2] & 0x3F;
            std::string lesName = getLesName(sat, lesId);
            int logicalChannelNo = inputFrame.decodedFrame[*pos + 3];
            int packetNo = inputFrame.decodedFrame[*pos + 4];
            //in bytes, presentation agnostic
            int j = *pos + 5;
            ret.payload.presentation = 0;
            for(int i = 0; j < *pos + ret.packetLength - 2; i++) {
                char c = inputFrame.decodedFrame[j];
                ret.payload.data8Bit.push_back(c);
                j++;
            }
            ret.payload.presentation = IsBinary(ret.payload.data8Bit, true) ? 7 : 0;
            ret.packetVars.insert(std::pair<std::string, std::string>("sat", std::to_string(sat)));
            ret.packetVars.insert(std::pair<std::string, std::string>("satName", satName));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesId", std::to_string(lesId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("lesName", lesName));
            ret.packetVars.insert(std::pair<std::string, std::string>("logicalChannelNo", std::to_string(logicalChannelNo)));
            ret.packetVars.insert(std::pair<std::string, std::string>("packetNo", std::to_string(packetNo)));
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_AB(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_COMPLETE;
            int lesListLength = inputFrame.decodedFrame[*pos + 1];
            std::ostringstream os;
            for(int i = 0; i < 1; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << (uint16_t)inputFrame.decodedFrame[*pos + 2 + i];
            }
            std::string stationStartHex = os.str();
            int stationCount = inputFrame.decodedFrame[*pos + 3];
            std::string stations = getStations(inputFrame.decodedFrame, stationCount, *pos + 4);
            ret.packetVars.insert(std::pair<std::string, std::string>("lesListLength", std::to_string(lesListLength)));
            ret.packetVars.insert(std::pair<std::string, std::string>("stationStartHex", stationStartHex));
            ret.packetVars.insert(std::pair<std::string, std::string>("stations", stations));
            ret.packetVars.insert(std::pair<std::string, std::string>("stationCount", std::to_string(stationCount)));
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_AC(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            return ret; //not implemented yet
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_AD(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            return ret; //not implemented yet
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_B1(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_PARTIAL;
            int messageType = inputFrame.decodedFrame[*pos + 2];
            std::string serviceCodeAndAddressName = getServiceCodeAndAddressName(messageType);
            int continuation = (inputFrame.decodedFrame[*pos + 3] & 0x80) >> 7;
            int priority = (inputFrame.decodedFrame[*pos + 3] & 0x60) >> 5;
            std::string priorityText = getPriority(priority);
            bool isDistress = priority == 3;
            int repetition = inputFrame.decodedFrame[*pos + 3] & 0x1F;
            int messageId = inputFrame.decodedFrame[*pos + 4] << 8 | inputFrame.decodedFrame[*pos + 5];
            int packetNo = inputFrame.decodedFrame[*pos + 6];
            bool isNewPayload = packetNo == 1;
            ret.payload.presentation = inputFrame.decodedFrame[*pos + 7];
            ret.packetVars.insert(std::pair<std::string, std::string>("messageType", std::to_string(messageType)));
            ret.packetVars.insert(std::pair<std::string, std::string>("serviceCodeAndAddressName", serviceCodeAndAddressName));
            ret.packetVars.insert(std::pair<std::string, std::string>("continuation", std::to_string(continuation)));
            ret.packetVars.insert(std::pair<std::string, std::string>("priority", std::to_string(priority)));
            ret.packetVars.insert(std::pair<std::string, std::string>("priorityText", priorityText));
            ret.packetVars.insert(std::pair<std::string, std::string>("isDistress", std::to_string(isDistress)));
            ret.packetVars.insert(std::pair<std::string, std::string>("repetition", std::to_string(repetition)));
            ret.packetVars.insert(std::pair<std::string, std::string>("messageId", std::to_string(messageId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("packetNo", std::to_string(packetNo)));
            ret.packetVars.insert(std::pair<std::string, std::string>("isNewPayload", std::to_string(isNewPayload)));
            int addressLength = getAddressLength(messageType);
            //NAV/MET coordinator... area... TODO
            if(*pos + 8 + addressLength >= inputFrame.length) {
                return ret;
            }
            uint8_t address[addressLength];
            std::copy(&inputFrame.decodedFrame[*pos+8], &inputFrame.decodedFrame[*pos+8+addressLength], address);
            std::ostringstream os;
            for(int i = 0; i < addressLength - 1; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << (uint16_t)address[i + 1];
            }
            std::string addressHex = os.str();
            int payloadLength = ret.packetLength - 2 - 8 - addressLength;
            int k = *pos + 8 + addressLength;
            for(int i = 0; k < *pos + 8 + addressLength + payloadLength; i++) {
                ret.payload.data8Bit.push_back(inputFrame.decodedFrame[k]);
                k++;
            }
            ret.packetVars.insert(std::pair<std::string, std::string>("addressHex", addressHex));
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_B2(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_PARTIAL;
            int messageType = inputFrame.decodedFrame[*pos + 2];
            std::string serviceCodeAndAddressName = getServiceCodeAndAddressName(messageType);
            int continuation = (inputFrame.decodedFrame[*pos + 3] & 0x80) >> 7;
            int priority = (inputFrame.decodedFrame[*pos + 3] & 0x60) >> 5;
            std::string priorityText = getPriority(priority);
            bool isDistress = priority == 3;
            int repetition = inputFrame.decodedFrame[*pos + 3] & 0x1F;
            int messageId = inputFrame.decodedFrame[*pos + 4] << 8 | inputFrame.decodedFrame[*pos + 5];
            int packetNo = inputFrame.decodedFrame[*pos + 6];
            bool isNewPayload = packetNo == 1;
            ret.payload.presentation = inputFrame.decodedFrame[*pos + 7];
            ret.packetVars.insert(std::pair<std::string, std::string>("messageType", std::to_string(messageType)));
            ret.packetVars.insert(std::pair<std::string, std::string>("serviceCodeAndAddressName", serviceCodeAndAddressName));
            ret.packetVars.insert(std::pair<std::string, std::string>("continuation", std::to_string(continuation)));
            ret.packetVars.insert(std::pair<std::string, std::string>("priority", std::to_string(priority)));
            ret.packetVars.insert(std::pair<std::string, std::string>("priorityText", priorityText));
            ret.packetVars.insert(std::pair<std::string, std::string>("isDistress", std::to_string(isDistress)));
            ret.packetVars.insert(std::pair<std::string, std::string>("repetition", std::to_string(repetition)));
            ret.packetVars.insert(std::pair<std::string, std::string>("messageId", std::to_string(messageId)));
            ret.packetVars.insert(std::pair<std::string, std::string>("packetNo", std::to_string(packetNo)));
            ret.packetVars.insert(std::pair<std::string, std::string>("isNewPayload", std::to_string(isNewPayload)));
            int addressLength = getAddressLength(messageType);
            //NAV/MET coordinator... area... TODO
            if(*pos + 8 + addressLength >= inputFrame.length) {
                return ret;
            }
            uint8_t address[addressLength];
            std::copy(&inputFrame.decodedFrame[*pos+8], &inputFrame.decodedFrame[*pos+8+addressLength], address);
            std::ostringstream os;
            for(int i = 0; i < addressLength - 1; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << (uint16_t)address[i + 1];
            }
            std::string addressHex = os.str();
            int payloadLength = ret.packetLength - 2 - 8 - addressLength;
            int k = *pos + 8 + addressLength;
            for(int i = 0; k < *pos + 8 + addressLength + payloadLength; i++) {
                ret.payload.data8Bit.push_back(inputFrame.decodedFrame[k]);
                k++;
            }
            ret.packetVars.insert(std::pair<std::string, std::string>("addressHex", addressHex));
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_BD(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_COMPLETE;
            int multiframePacketLength = 0;
            int multiframePacketDescriptor = inputFrame.decodedFrame[*pos + 2] & 0xFF;
            std::ostringstream os;
            os << std::setfill('0') << std::setw(2) << std::right << std::hex << multiframePacketDescriptor;
            std::string multiframePacketDescriptorHex = os.str();
            if(multiframePacketDescriptor >> 7 == 0) {
                // Short packet descriptor
                // The packet length including CRC does not include byte 0, we add 1
                multiframePacketLength = (multiframePacketDescriptor & 0x0F) + 1;
            } else if(multiframePacketDescriptor >> 6 == 0x02) {
                /// Medium packet descriptor
                /// The packet length including CRC does not include the first 2 bytes, we add 2
                multiframePacketLength = inputFrame.decodedFrame[*pos + 3] + 2;
            }
            //compose encapsulated packet, length includes the CRC for the composed packet
            ret.mfp.packetData.resize(multiframePacketLength);
            //set the length of payload data to be used for the new packet
            //the payload starts at index + 2 and does not include the CRC
            ret.mfp.firstPartCount = ret.packetLength - 2 - 2;
            ret.mfp.isMFP = true;
            std::copy(&inputFrame.decodedFrame[*pos + 2], &inputFrame.decodedFrame[*pos + 2 + ret.mfp.firstPartCount], ret.mfp.packetData.begin());
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDecoder::decode_BE(decoder::Decoder::decoder_result inputFrame, int* pos, PacketDecoder::packetDecoder_multiFramePacket* mfa) {
            PacketDecoder::packetDecoder_result ret = basicDecode(inputFrame, pos);
            if(!ret.isCrc) {
                return ret;
            }
            ret.decodingStage = PACKETDECODER_DECODING_STAGE_COMPLETE;
            //2 for CRC
            //2 for starting of packet
            int j = *pos + 2;
            int actualLength = 0;
            for(int i = 0; j < *pos + ret.packetLength - 2; i++) {
                ret.payload.data8Bit.push_back(inputFrame.decodedFrame[j]);
                if((int)mfa->packetData.size() > i + mfa->firstPartCount) {
                    mfa->packetData[i + mfa->firstPartCount] = (inputFrame.decodedFrame[j]);
                }
                j++;
                actualLength++;
            }
            mfa->firstPartCount += actualLength;
            mfa->isReady = mfa->firstPartCount >= (int)mfa->packetData.size() - 2;
            return ret;
        }
        int PacketDecoder::computeCRC(uint8_t *decodedFrame, int pos, int length) {
            short C0 = 0;
            short C1 = 0;
            uint8_t CB1;
            uint8_t CB2;
            uint8_t B;
            int i = 0;
            while (i < length) {
                if (i < length - 2) {
                    B = decodedFrame[pos + i];
                } else {
                    B = 0;
                }
                C0 += B;
                C1 += C0;
                i++;
            }
            CB1 = (uint8_t)(C0 - C1);
            CB2 = (uint8_t)(C1 - 2 * C0);
            return (CB1 << 8) | CB2;
        }
        std::string PacketDecoder::getSatName(int sat) {
            switch (sat) {
                case 0:
                    return "Atlantic Ocean Region West (AOR-W)";
                case 1:
                    return "Atlantic Ocean Region East (AOR-E)";
                case 2:
                    return "Pacific Ocean Region (POR)";
                case 3:
                    return "Indian Ocean Region (IOR)";
                case 9:
                    return "All Ocean Regions Covered by the LES";
                default:
                    return "Unknown";
            }
        }
        std::string PacketDecoder::getLesName(int sat, int lesId) {
            int value = lesId + sat * 100;
            std::string name;
            switch (value) {
                case 001:
                case 101:
                case 201:
                case 301:
                    name = "Vizada-Telenor, USA";
                    break;

                case 002:
                case 102:
                case 302:
                    name = "Stratos Global (Burum-2), Netherlands";
                    break;

                case 202:
                    name = "Stratos Global (Aukland), New Zealand";
                    break;

                case 003:
                case 103:
                case 203:
                case 303:
                    name = "KDDI Japan";
                    break;

                case 004:
                case 104:
                case 204:
                case 304:
                    name = "Vizada-Telenor, Norway";
                    break;

                case 044:
                case 144:
                case 244:
                case 344:
                    name = "NCS";
                    break;

                case 105:
                case 335:
                    name = "Telecom, Italia";
                    break;

                case 305:
                case 120:
                    name = "OTESTAT, Greece";
                    break;

                case 306:
                    name = "VSNL, India";
                    break;

                case 110:
                case 310:
                    name = "Turk Telecom, Turkey";
                    break;

                case 211:
                case 311:
                    name = "Beijing MCN, China";
                    break;

                case 012:
                case 112:
                case 212:
                case 312:
                    name = "Stratos Global (Burum), Netherlands";
                    break;

                case 114:
                    name = "Embratel, Brazil";
                    break;

                case 116:
                case 316:
                    name = "Telekomunikacja Polska, Poland";
                    break;

                case 117:
                case 217:
                case 317:
                    name = "Morsviazsputnik, Russia";
                    break;

                case 021:
                case 121:
                case 221:
                case 321:
                    name = "Vizada (FT), France";
                    break;

                case 127:
                case 327:
                    name = "Bezeq, Israel";
                    break;

                case 210:
                case 328:
                    name = "Singapore Telecom, Singapore";
                    break;

                case 330:
                    name = "VISHIPEL, Vietnam";
                    break;

                default:
                    name = "Unknown";
                    break;
            }
            return std::to_string(value) + ", " + name;
        }
        std::string PacketDecoder::getServiceCodeAndAddressName(int code) {
            switch (code) {
                case 0x00:
                    return "System, All ships (general call)";
                case 0x02:
                    return "FleetNET, Group Call";
                case 0x04:
                    return "SafetyNET, Navigational, Meteorological or Piracy Warning to a Rectangular Area";
                case 0x11:
                    return "System, Inmarsat System Message";
                case 0x13:
                    return "SafetyNET, Navigational, Meteorological or Piracy Coastal Warning";
                case 0x14:
                    return "SafetyNET, Shore-to-Ship Distress Alert to Circular Area";
                case 0x23:
                    return "System, EGC System Message";
                case 0x24:
                    return "SafetyNET, Navigational, Meteorological or Piracy Warning to a Circular Area";
                case 0x31:
                    return "SafetyNET, NAVAREA/METAREA Warning, MET Forecast or Piracy Warning to NAVAREA/METAREA";
                case 0x33:
                    return "System, Download Group Identity";
                case 0x34:
                    return "SafetyNET, SAR Coordination to a Rectangular Area";
                case 0x44:
                    return "SafetyNET, SAR Coordination to a Circular Area";
                case 0x72:
                    return "FleetNET, Chart Correction Service";
                case 0x73:
                    return "SafetyNET, Chart Correction Service for Fixed Areas";
                default:
                    return "Unknown";
            }
        }
        std::string PacketDecoder::getPriority(int priority) {
            switch (priority) {
                case -1:
                    return "Message";
                case 0:
                    return "Routine";
                case 1:
                    return "Safety";
                case 2:
                    return "Urgency";
                case 3:
                    return "Distress";
                default:
                    return "Unknown";
            }
        }
        int PacketDecoder::getAddressLength(int messageType) {
            switch (messageType) {
                case 0x00:
                    return 3;
                case 0x11:
                case 0x31:
                    return 4;
                case 0x02:
                case 0x72:
                    return 5;
                case 0x13:
                case 0x23:
                case 0x33:
                case 0x73:
                    return 6;
                case 0x04:
                case 0x14:
                case 0x24:
                case 0x34:
                case 0x44:
                    return 7;
                default:
                    return 3;
            }
        }
        bool PacketDecoder::IsBinary(std::vector<uint8_t> data, bool checkAll) {
            bool isBinary = false;
            //try first 13 characters if not check all
            int check = 13;
            if (!checkAll) {
                check = std::min(check, (int)data.size()-2);
            } else {
                check = data.size();
            }
            for (int i = 0; i < check; i++) {
                if (!isBinary) {
                    char chr = (char)(data[i] & 0x7F);
                    switch (chr) {
                        case 0x01:
                        case 0x05:
                        case 0x06:
                        case 0x07:
                        case 0x08:
                        case 0x0B:
                        case 0x0C:
                        case 0x0E:
                        case 0x0F:
                        case 0x10:
                        case 0x11:
                        case 0x12:
                        case 0x13:
                        case 0x14:
                        case 0x15:
                        case 0x16:
                        case 0x17:
                        case 0x18:
                        case 0x19:
                        case 0x1A:
                        case 0x1C:
                        case 0x1D:
                        case 0x1E:
                        case 0x1F:
                        //case 0xA4:
                        case '$':
                            isBinary = true;
                            break;
                    }
                    if (isBinary)
                    {
                        break;
                    }
                }
            }
            return isBinary;
        }
        std::string PacketDecoder::getStations(uint8_t data[], int stationCount, int pos) {
            std::string stations;
            int j = pos;
            for(int i=0; i < stationCount; i++) {
                stations += "Station: " + std::to_string(i) + "\n";
                int sat = data[j] >> 6 & 0x03;
                stations += "   sat: " + std::to_string(sat) + "\n";
                std::string satName = getSatName(sat);
                stations += "   satName: " + satName + "\n";
                int lesId = data[j] & 0x3F;
                stations += "   lesId: " + std::to_string(lesId) + "\n";
                std::string lesName = getLesName(sat, lesId);
                stations += "   lesName: " + lesName + "\n";
                int servicesStart = data[j + 1];
                stations += "   servicesStart: " + std::to_string(servicesStart) + "\n";
                int iss = data[j + 2] << 8 | data[j + 3];
                stations += "   Services:";
                stations += "       MaritimeDistressAlerting: " + std::to_string((iss & 0x8000) >> 15 == 1) + "\n";
                stations += "       SafetyNet: " + std::to_string((iss & 0x4000) >> 14 == 1) + "\n";
                stations += "       InmarsatC: " + std::to_string((iss & 0x2000) >> 13 == 1) + "\n";
                stations += "       StoreFwd: " + std::to_string((iss & 0x1000) >> 12 == 1) + "\n";
                stations += "       HalfDuplex: " + std::to_string((iss & 0x800) >> 11 == 1) + "\n";
                stations += "       FullDuplex: " + std::to_string((iss & 0x400) >> 10 == 1) + "\n";
                stations += "       ClosedNetwork: " + std::to_string((iss & 0x200) >> 9 == 1) + "\n";
                stations += "       FleetNet: " + std::to_string((iss & 0x100) >> 8 == 1) + "\n";
                stations += "       PrefixSF: " + std::to_string((iss & 0x80) >> 7 == 1) + "\n";
                stations += "       LandMobileAlerting: " + std::to_string((iss & 0x40) >> 6 == 1) + "\n";
                stations += "       AeroC: " + std::to_string((iss & 0x20) >> 5 == 1) + "\n";
                stations += "       ITA2: " + std::to_string((iss & 0x10) >> 4 == 1) + "\n";
                stations += "       DATA: " + std::to_string((iss & 0x08) >> 3 == 1) + "\n";
                stations += "       BasicX400: " + std::to_string((iss & 0x04) >> 2 == 1) + "\n";
                stations += "       EnhancedX400: " + std::to_string((iss & 0x02) >> 1 == 1) + "\n";
                stations += "       LowPowerCMES: " + std::to_string((iss & 0x01) == 1);
                double downlinkChannelMhz = ((data[j + 4] << 8 | data[j + 4]) - 8000) * 0.0025 + 1530.5;
                stations += "  downlinkChannelMhz : " + std::to_string(downlinkChannelMhz) + "\n";
                j+= 6;
            }
            return stations;
        }
        std::string PacketDecoder::getServices_short(uint8_t is8) {
            std::string services;
            services += "MaritimeDistressAlerting: " + std::to_string((is8 & 0x80) >> 7 == 1) + "\n";
            services += "SafetyNet: " + std::to_string((is8 & 0x40) >> 6 == 1) + "\n";
            services += "InmarsatC: " + std::to_string((is8 & 0x20) >> 5 == 1) + "\n";
            services += "StoreFwd: " + std::to_string((is8 & 0x10) >> 4 == 1) + "\n";
            services += "HalfDuplex: " + std::to_string((is8 & 8) >> 3 == 1) + "\n";
            services += "FullDuplex: " + std::to_string((is8 & 4) >> 2 == 1) + "\n";
            services += "ClosedNetwork: " + std::to_string((is8 & 2) >> 1 == 1) + "\n";
            services += "FleetNet: " + std::to_string((is8 & 1) == 1);
            return services;
        }
        std::string PacketDecoder::getServices(int iss) {
            std::string services;
            services += "MaritimeDistressAlerting: " + std::to_string((iss & 0x8000) >> 15 == 1) + "\n";
            services += "SafetyNet: " + std::to_string((iss & 0x4000) >> 14 == 1) + "\n";
            services += "InmarsatC: " + std::to_string((iss & 0x2000) >> 13 == 1) + "\n";
            services += "StoreFwd: " + std::to_string((iss & 0x1000) >> 12 == 1) + "\n";
            services += "HalfDuplex: " + std::to_string((iss & 0x800) >> 11 == 1) + "\n";
            services += "FullDuplex: " + std::to_string((iss & 0x400) >> 10 == 1) + "\n";
            services += "ClosedNetwork: " + std::to_string((iss & 0x200) >> 9 == 1) + "\n";
            services += "FleetNet: " + std::to_string((iss & 0x100) >> 8 == 1) + "\n";
            services += "PrefixSF: " + std::to_string((iss & 0x80) >> 7 == 1) + "\n";
            services += "LandMobileAlerting: " + std::to_string((iss & 0x40) >> 6 == 1) + "\n";
            services += "AeroC: " + std::to_string((iss & 0x20) >> 5 == 1) + "\n";
            services += "ITA2: " + std::to_string((iss & 0x10) >> 4 == 1) + "\n";
            services += "DATA: " + std::to_string((iss & 0x08) >> 3 == 1) + "\n";
            services += "BasicX400: " + std::to_string((iss & 0x04) >> 2 == 1) + "\n";
            services += "EnhancedX400: " + std::to_string((iss & 0x02) >> 1 == 1) + "\n";
            services += "LowPowerCMES: " + std::to_string((iss & 0x01) == 1);
            return services;
        }
        std::string PacketDecoder::getDescriptorAsText(uint8_t descriptor_b) {
            std::string descriptor;
            switch(descriptor_b) {
                case 0x27:
                    descriptor = "Logical Channel Clear";
                    break;
                case 0x2A:
                    descriptor = "Inbound Message Ack";
                    break;
                case 0x08:
                    descriptor = "Acknowledgement Request";
                    break;
                case 0x6C:
                    descriptor = "Signalling Channel";
                    break;
                case 0x7D:
                    descriptor = "Bulletin Board";
                    break;
                case 0x81:
                    descriptor = "Announcement";
                    break;
                case 0x83:
                    descriptor = "Logical Channel Assignment";
                    break;
                case 0x91:
                    descriptor = "Distress Alert Ack";
                    break;
                case 0x92:
                    descriptor = "Login Ack";
                    break;
                case 0x9A:
                    descriptor = "Enhanced Data Report Ack";
                    break;
                case 0xA0:
                    descriptor = "Distress Test Request";
                    break;
                case 0xA3:
                    descriptor = "Individual Poll";
                    break;
                case 0xA8:
                    descriptor = "Confirmation";
                    break;
                case 0xAA:
                    descriptor = "Message";
                    break;
                case 0xAB:
                    descriptor = "Les List";
                    break;
                case 0xAC:
                    descriptor = "Request Status";
                    break;
                case 0xAD:
                    descriptor = "Test Result";
                    break;
                case 0xB1:
                    descriptor = "EGC double header, part 1";
                    break;
                case 0xB2:
                    descriptor = "EGC double header, part 2";
                    break;
                case 0xBD:
                    descriptor = "Multiframe Packet Start";
                    break;
                case 0xBE:
                    descriptor = "Multiframe Packet Continue";
                    break;
                default:
                    descriptor = "Unknown";
                    break;
            }
            return descriptor;
        }
        //END CLASS PacketDecoder

        //START CLASS PacketDetector
        PacketDetector::PacketDetector() {
            packetDecoder = new PacketDecoder();
        }
        std::vector<PacketDecoder::packetDecoder_result> PacketDetector::process(decoder::Decoder::decoder_result inputFrame) {
            std::vector<PacketDecoder::packetDecoder_result> ret;
            int pos = 0;
            do {
                PacketDecoder::packetDecoder_result res = detect(inputFrame, &pos);
                ret.push_back(res);
            } while(pos > 0);
            return ret;
        }
        PacketDecoder::packetDecoder_result PacketDetector::detect(decoder::Decoder::decoder_result inputFrame, int* pos) {
            PacketDecoder::packetDecoder_result ret;
            ret.isDecodedPacket = false;
            if(*pos >= inputFrame.length) {
                *pos = 0;
                return ret;
            }
            switch(inputFrame.decodedFrame[*pos]) {
                //No more data
                case 0x00:
                    *pos = 0;
                    return ret;
                //27 - Logical Channel Clear
                case 0x27:
                    ret = packetDecoder->decode_27(inputFrame, pos);
                    break;
                //2A - Inbound Message Ack.
                case 0x2A:
                    ret = packetDecoder->decode_2A(inputFrame, pos);
                    break;
                //08 - Acknowledgement Request
                case 0x08:
                    ret = packetDecoder->decode_08(inputFrame, pos);
                    break;
                //6C - Signalling Channel
                case 0x6C:
                    ret = packetDecoder->decode_6C(inputFrame, pos);
                    break;
                //7D - Bulletin Board
                case 0x7D:
                    ret = packetDecoder->decode_7D(inputFrame, pos);
                    break;
                //81 - Announcement
                case 0x81:
                    ret = packetDecoder->decode_81(inputFrame, pos);
                    break;
                //83 - Logical Channel Assignment
                case 0x83:
                    ret = packetDecoder->decode_83(inputFrame, pos);
                    break;
                //91 - Distress Alert Ack.
                case 0x91:
                    ret = packetDecoder->decode_91(inputFrame, pos);
                    break;
                //92 - Login Ack.
                case 0x92:
                    ret = packetDecoder->decode_92(inputFrame, pos);
                    break;
                //9A - Enhanced Data Report Ack.
                case 0x9A:
                    ret = packetDecoder->decode_9A(inputFrame, pos);
                    break;
                //A0 - Distress Test Request
                case 0xA0:
                    ret = packetDecoder->decode_A0(inputFrame, pos);
                    break;
                //A3 - Individual Poll
                case 0xA3:
                    ret = packetDecoder->decode_A3(inputFrame, pos);
                    break;
                //A8 - Confirmation
                case 0xA8:
                    ret = packetDecoder->decode_A8(inputFrame, pos);
                    break;
                //AA - Message
                case 0xAA:
                    ret = packetDecoder->decode_AA(inputFrame, pos);
                    break;
                //AB - Les List
                case 0xAB:
                    ret = packetDecoder->decode_AB(inputFrame, pos);
                    break;
                //AC - Request Status
                case 0xAC:
                    ret = packetDecoder->decode_AC(inputFrame, pos);
                    break;
                //AD - Test Result
                case 0xAD:
                    ret = packetDecoder->decode_AD(inputFrame, pos);
                    break;
                //B1 - EGC double header, part 1 (16 bytes of data)
                case 0xB1:
                    ret = packetDecoder->decode_B1(inputFrame, pos);
                    break;
                //B2 - EGC double header, part 2
                case 0xB2:
                    ret = packetDecoder->decode_B2(inputFrame, pos);
                    break;
                //BD - Multiframe Packet
                case 0xBD:
                    ret = packetDecoder->decode_BD(inputFrame, pos);
                    multiStreamFrameElements = ret.mfp;
                    break;
                //BE - Multiframe Packet Continue
                case 0xBE:
                    ret = packetDecoder->decode_BE(inputFrame, pos, &multiStreamFrameElements);
                    if(!multiStreamFrameElements.isMFP) {
                        //if BD is not received before, don't try to decode BE
                        break;
                    }
                    if(multiStreamFrameElements.isReady) {
                        //decode
                        decoder::Decoder::decoder_result ndfa;
                        for(int i = 0; i < (int)multiStreamFrameElements.packetData.size(); i++) {
                            ndfa.decodedFrame[i] = multiStreamFrameElements.packetData[i];
                        }
                        ndfa.length = multiStreamFrameElements.packetData.size();
                        ndfa.frameNumber = inputFrame.frameNumber;
                        ndfa.BER = inputFrame.BER;
                        ndfa.isHardDecision = inputFrame.isHardDecision;
                        ndfa.isMidStreamReversePolarity = inputFrame.isMidStreamReversePolarity;
                        ndfa.isReversedPolarity = inputFrame.isReversedPolarity;
                        ndfa.isUncertain = inputFrame.isUncertain;
                        ndfa.timestamp = inputFrame.timestamp;
                        //recursively decode
                        int pos1 = 0;
                        ret = detect(ndfa, &pos1);
                        multiStreamFrameElements.isMFP = false;
                    }
                    break;
                default:
                    ret = packetDecoder->basicDecode(inputFrame, pos);
                    break;
            }
            *pos += ret.packetLength;
            ret.isDecodedPacket = true;
            return ret;
        }
        //END CLASS PacketDetector

        //START CLASS FrameParser
        FrameParser::FrameParser() {
            packetDetector = new PacketDetector();
        }
        std::vector<FrameParser::frameParser_result> FrameParser::parseFrame(decoder::Decoder::decoder_result inputFrame) {
            std::vector<frameParser_result> ret;
            std::vector<PacketDecoder::packetDecoder_result> detRes = packetDetector->process(inputFrame);
            for(int i = 0; i < (int)detRes.size();i++) {
                if(detRes[i].isDecodedPacket) {
                    frameParser_result res;
                    res.decoding_result = detRes[i];
                    ret.push_back(res);
                }
            }
            return ret;
        }
        //END CLASS FrameParser
    }
}
