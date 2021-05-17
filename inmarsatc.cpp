#include "inmarsatc.h"

namespace inmarsatc {
    namespace demodulator {
        //START CLASS FIR
        FIR::FIR(double b[], int len) {
            this->b = new double[len];
            std::copy(b, &b[len], this->b);
            psLength = len + 1;
            prevSamples = new double[len + 1];
        }
        double FIR::filter(double sample) {
            //std::rotate(prevSamples, (prevSamples.begin() + 1), prevSamples.end());
            for (int i = 0; i < psLength-1; i++) {
                prevSamples[i] = prevSamples[i + 1];
            }
            prevSamples[psLength - 1] = sample;
            return filter();
        }
        double FIR::filter() {
            int m = psLength - 1;
            int n = psLength;
            double y[n];
            for(int yi = 0; yi < n; yi++) {
                double t = 0.0;
                for(int bi = (m - 1); bi >= 0; bi--) {
                    if(yi - bi < 0) {
                        continue;
                    }
                    t += b[bi] * prevSamples[yi - bi];
                }
                y[yi] = t;
            }
            return y[psLength - 1];
        }
        void FIR::resetFilter() {
            prevSamples[0] = 0;
            prevSamples[1] = 0;
            prevSamples[2] = 0;
        }
        //END CLASS FIR

        //START CLASS RRC
        RRC::RRC(double alpha, int firstSize, double sampleRate, double symbolRate) {
            if((firstSize % 2) == 0) {
                firstSize++;
            }
            points = new double[firstSize];
            double t = sampleRate / symbolRate;
            double fi;
            for(int i = 0; i < firstSize; i++) {
                if(i == ((firstSize - 1) / 2)) {
                    points[i] = (4.0 * alpha + M_PI - M_PI * alpha) / (M_PI * sqrt(t));
                } else {
                    fi = (((double)i) - ((double)(firstSize - 1)) / 2.0);
                    if(abs(1.0 - pow(4.0 * alpha * fi / t, 2)) < 0.0000000001) {
                        points[i] = (alpha * ((M_PI - 2.0) * cos(M_PI / (4.0 * alpha)) + (M_PI + 2.0) * sin(M_PI / (4.0 * alpha))) / (M_PI * sqrt(2.0 * t)));
                    } else {
                        points[i] = (4.0 * alpha / (M_PI * sqrt(t)) * (cos((1.0 + alpha) * M_PI * fi / t) + t / (4.0 * alpha * fi) * sin((1.0 - alpha) * M_PI * fi / t)) / (1.0 - pow(4.0 * alpha * fi / t, 2)));
                    }
                }
            }
            rrcFIR = new FIR(points, firstSize);
        }
        double RRC::filter(double sample) {
            return rrcFIR->filter(sample);
        }
        //END CLASS RRC

        //START CLAS Gardner
        Gardner::Gardner(double sampleRate, double symbolRate) {
            ts = sampleRate / symbolRate;
            symbol2xTimer = 0;
            samplesAgo = 0;
            maxErrorxAggression = ts / 4;
            isOnPoint = 1;
            sampleLastOn = std::complex<double>(0, 0);
            aggression = 0.1;
            dummyOutput = std::complex<double>(0, 0);
        }
        bool Gardner::step(double re, double im, std::complex<double>* output, bool* onPoint) {
            //default return and output
            result = false;
            *output = dummyOutput;
            *onPoint = false;
            // these are samples after the carrier tracking and the RRC filter
            // current_sample = data_re(i) + 1i* data_im(i);
            currentSample = std::complex<double>(re, im);
            // sample at 2x bit rate
            symbol2xTimer += 2;
            // simple way to avoid double sampling due to timing jitter
            samplesAgo++;
            if(symbol2xTimer >= ts && complexMath::biggerThan(samplesAgo, (maxErrorxAggression + std::complex<double>(1, 0)))) {
                symbol2xTimer -= ts;
                samplesAgo = 0;
                isOnPoint = 1 - isOnPoint;
                if(isOnPoint == 1) {
                    // on time here
                    sampleThisOn = currentSample;
                    // calculate carrier timing error
                    error = (sampleThisOn - sampleLastOn) * sampleOff;
                    // steer symbol timing
                    symbol2xTimer = symbol2xTimer + (error * aggression).real();
                    // save on symbol
                    result = true;
                    *output = sampleThisOn;
                    *onPoint = true;
                    // save this on time for next on time
                    sampleLastOn = sampleThisOn;
                } else {
                    // we are off time here
                    sampleOff = currentSample;
                    // save off symbol
                    // we can control here how ofter we want the CMA to update its filter coefs
                    // if we set it to true, CMA will go for every sample not on point
                    result = true;
                    *output = sampleOff;
                    *onPoint = false;
                }
            }
            return result;
        }
        //END CLASS Gardner

        //START CLASS CMA
        CMA::CMA() {
            dummyOutput = std::complex<double>(0,0);
            // no of points of the fir filter
            for (int i = 0; i < CMA_SZ; i++) {
                cmaW[i] = std::complex<double>(0, 0);
            }
            //cmaW((cmaSz+1)/2)=1.0+1i*0; verify
            cmaW[(CMA_SZ) / 2] = std::complex<double>(1, 0);
            for (int i = 0; i < CMA_SZ; i++) {
                cmaX[i] = std::complex<double>(0, 0);
            }
            // there is some stats with this number and depends on what agc is set to
            beta = sqrt(2.0);
            //beta = 2;
            // this makes a huge difference
            //stepsize = 0.001;
            stepSize = 0.001;
        }
        bool CMA::step(std::complex<double> sample, bool isOnPoint, std::complex<double> *output) {
            //default return and output
            *output = dummyOutput;
            result = false;
            // load another sample these samples are running at 2x speed
            // but could be running faster if you want and i think would produce a
            // better equalizer but uses more cpu
            if(isOnPoint) {
                //AGC
                mean = (abs(sample) + abs(*output)) / 2.0;
                sample = sample / std::complex<double>(mean, 0);
            }
            // run samples though the equalizer
            // how often this runs is not critical but it should be constant wrt
            // symbol timing
            for(int i = 0; i < CMA_SZ - 1; i++) {
                cmaX[i] = cmaX[i + 1];
            }
            cmaX[CMA_SZ - 1] = sample;
            cmaEqualizerOut = std::complex<double>(0, 0);
            for(int i = 0; i < CMA_SZ; i++) {
                cmaEqualizerOut = cmaEqualizerOut + (cmaW[i] * cmaX[i]);
            }
            // this must be run only when you think cma_equalizer_out is an on point
            // if on time then update cma algo equalizer coeffs
            if(isOnPoint) {
                // save samples that come out of the equalizer that are on point
                *output = cmaEqualizerOut;
                result = true;
                // calc directional error. this one here only does modulus not
                // rotation.the rotation one is in jdsca but is still simple
                // error = cma_equalizer_out * ((abs(cma_equalizer_out)) ^ 2 - beta); verify
                error = cmaEqualizerOut * std::complex<double>(pow(abs(cmaEqualizerOut), 2) - beta, 0.0);
                if(abs(error) == 0) {
                    error.real(1); //avoid stuck if samples are equal to zero
                }
                // step filer coeffs in direction opposite direction of error
                for(int i = 0; i < CMA_SZ; i++) {
                    cmaW[i] -= stepSize * error * std::conj(cmaX[i]);
                }
            }
            if(!std::isfinite(abs(*output))) {
                *output=dummyOutput;
                CMAReset();
            }
            return result;
        }
        void CMA::CMAReset() {
            for(int i = 0; i < CMA_SZ; i++) {
                cmaW[i] = std::complex<double>(0, 0);
            }
            //cmaW((cmaSz+1)/2)=1.0+1i*0; verify
            cmaW[(CMA_SZ) / 2] = std::complex<double>(1, 0);
            for(int i = 0; i < CMA_SZ; i++) {
                cmaX[i] = std::complex<double>(0, 0);
            }
        }
        //END CLASS CMA

        //START CLASS AGC
        AGC::AGC() {
            mAGC = 0.0;
        }
        void AGC::apply(double* i, double* q) {
            double magnitude = sqrt((*i * *i) + (*q * *q));
            if(magnitude > mAGC) {
                mAGC = (1.0 - 1.0 / 250.0) * mAGC + (1.0 / 250.0) * magnitude;
            } else {
                mAGC = (1.0 - 1.0 / 1000.0) * mAGC + (1.0 / 1000.0) * magnitude;
            }
            if(mAGC >= 1.0) {
                *i = *i / mAGC;
                *q = *q / mAGC;
            }
        }
        void AGC::apply(std::complex<double> *value) {
            double magnitude = sqrt((value->real() * value->real()) + (value->imag() * value->imag()));
            if (magnitude > mAGC) {
                mAGC = (1.0 - 1.0 / 250.0) * mAGC + (1.0 / 250.0) * magnitude;
            } else {
                mAGC = (1.0 - 1.0 / 1000.0) * mAGC + (1.0 / 1000.0) * magnitude;
            }
            if (mAGC >= 1.0) {
                *value = std::complex<double>(value->real() / mAGC, value->imag() / mAGC);
            }
        }
        void AGC::reset() {
            mAGC = 0.0;
        }
        //END CLASS AGC

        //START CLASS Demodulator
        Demodulator::Demodulator() {
            flagcounter = 0;
            // carrier
            omega = 2.0 * M_PI * freq / DEMODULATOR_SAMPLERATE;
            // I and Q LPFs
            lpf1 = new FIR(blpf12, 2);
            lpf2 = new FIR(blpf12, 2);
            // loop LPF
            lpf3 = new FIR(blpf3, 2);
            rrcRe = new RRC(1, 17, DEMODULATOR_SAMPLERATE, DEMODULATOR_SYMBOLRATE);
            rrcIm = new RRC(1, 17, DEMODULATOR_SAMPLERATE, DEMODULATOR_SYMBOLRATE);
            rrcRe2 = new RRC(1, 11, DEMODULATOR_SAMPLERATE, DEMODULATOR_SYMBOLRATE);
            rrcIm2 = new RRC(1, 11, DEMODULATOR_SAMPLERATE, DEMODULATOR_SYMBOLRATE);
            gardner = new Gardner(DEMODULATOR_SAMPLERATE, DEMODULATOR_SYMBOLRATE);
            cma = new CMA();
            cmaEnabled = true;
            agc = new AGC();
            agcEnabled = true;
            loFreq = DEMODULATOR_LOW_FREQ;
            hiFreq = DEMODULATOR_HIGH_FREQ;
        }
        bool Demodulator::isCmaEnabled() {
            return cmaEnabled;
        }
        bool Demodulator::isAgcEnabled() {
            return agcEnabled;
        }
        int Demodulator::getLowFreq() {
            return loFreq;
        }
        int Demodulator::getHighFreq() {
            return hiFreq;
        }
        double Demodulator::getCenterFreq() {
            return freq;
        }
        bool Demodulator::getIsInSync() {
            return isInSync;
        }
        int Demodulator::getNoSyncCount() {
            return noSyncCounter;
        }
        std::complex<double> Demodulator::getScatterPoint() {
            return lastScatterPoint;
        }
        void Demodulator::setCmaEnabled(bool cmaEnabled) {
            this->cmaEnabled = cmaEnabled;
        }
        void Demodulator::setAgcEnabled(bool agcEnabled) {
            this->agcEnabled = agcEnabled;
        }
        void Demodulator::setLowFreq(int lowFreq) {
            this->loFreq = lowFreq;
        }
        void Demodulator::setHighFreq(int highFreq) {
            this->hiFreq = highFreq;
        }
        void Demodulator::setCenterFreq(double centerFreq) {
            freq_mtx.try_lock();
            this->freq = centerFreq;
            freq_mtx.unlock();
        }
        void Demodulator::cmaReset() {
            cma->CMAReset();
        }
        std::vector<Demodulator::demodulator_result> Demodulator::demodulate(std::complex<double>* samples, int length) {
            std::vector<demodulator_result> ret;
            double vI;
            double vQ;
            double I;
            double Q;
            double magnitude = 0.0;
            double meanMagnitude;
            double Irrc;
            double Qrrc;
            double Irrc2;
            double Qrrc2;
            // sample normalization
            for(int i = 0; i < length; i++) {
                magnitude += abs(samples[i]);
            }
            meanMagnitude = magnitude / length;
            if(meanMagnitude == 0) {
                meanMagnitude = 1;
            }
            for(int i = 0; i < length; i++) {
                samples[i] = samples[i] / meanMagnitude;
            }
            //process current set of samples
            double err_summed = 0;
            for(int i = 0; i < length; i++) {
                freq_mtx.try_lock();
                // Costas classic carrier recovery, ignore/eliminate the phase info)
                // Use both sin and cos from the VCO and create the voltage by combining them
                freq = (omega * DEMODULATOR_SAMPLERATE) / (2.0 * M_PI);
                // phase
                phase = phase + omega + alpha * error;
                // carrier
                omega = omega + beta * error;
                // keep frequency in range
                if(freq < loFreq) {
                    freq = loFreq;
                    omega = 2.0 * M_PI * freq / DEMODULATOR_SAMPLERATE;
                }
                if(freq > hiFreq) {
                    freq = hiFreq;
                    omega = 2.0 * M_PI * freq / DEMODULATOR_SAMPLERATE;
                }
                freq_mtx.unlock();
                // keep phase in range
                while(phase > 2.0 * M_PI) {
                    phase -= 2.0 * M_PI;
                }
                // VCO
                vI = cos(phase);
                vQ = -sin(phase);
                // mixer
                I = vI * samples[i].real();
                Q = vQ * samples[i].imag();
                // LPFs
                I = lpf1->filter(I);
                Q = lpf2->filter(Q);
                // VCO error
                error = I * Q;
                // loop filter
                error = lpf3->filter(error);
                // Summing up all errors in this args sample set
                // As the carrier gets locked, the positive and the negative error are equally distributed and
                // their sum gets closer and closer to zero.
                err_summed += error;
                // RRC
                // We are trying to get to make every single pulse a sharp triangle to provide the input for the Gardner
                Irrc = rrcRe->filter(I);
                Qrrc = rrcIm->filter(Q);
                Irrc2 = rrcRe2->filter(Irrc);
                Qrrc2 = rrcIm2->filter(Qrrc);
                // Gardner
                bool isOnPoint;
                if(!gardner->step(Irrc2, Qrrc2, &gardnerOutputSample, &isOnPoint)) {
                    continue;
                }
                // CMA
                if(cmaEnabled) {
                    if(agcEnabled) {
                        agc->apply(&gardnerOutputSample);
                    }
                    if(!cma->step(gardnerOutputSample, isOnPoint, &cmaOutputSample)) {
                        continue;
                    }
                }
                if(!isOnPoint) {
                    continue;
                }
                if(cmaEnabled) {
                    scatterPoint = cmaOutputSample;
                } else {
                    scatterPoint = gardnerOutputSample;
                }
                //the I output is the demodulated soft symbol
                if(std::isnan(scatterPoint.real())) {
                    continue;
                }
                lastScatterPoint = scatterPoint;
                //Determine value sign
                uint8_t iBit = (((scatterPoint.real() > 0.0) - (scatterPoint.real() < 0.0)) + 1.0) / 2.0;
                //rotate right
                std::copy(symbolBuffer, &symbolBuffer[DEMODULATOR_SYMBOLSPERCHUNK], &symbolBuffer[1]);
                symbolBuffer[0] = iBit;
                flagcounter++;
                if(flagcounter == DEMODULATOR_SYMBOLSPERCHUNK) {
                    demodulator_result res;
                    std::copy(symbolBuffer, &symbolBuffer[DEMODULATOR_SYMBOLSPERCHUNK], res.bitsDemodulated);
                    std::reverse(res.bitsDemodulated, res.bitsDemodulated + DEMODULATOR_SYMBOLSPERCHUNK);
                    res.meanMagnitude = meanMagnitude;
                    ret.push_back(res);
                    flagcounter = 0;
                }
            }
            isInSync = (std::abs(err_summed) < DEMODULATOR_MAXABSERRSUM);
            noSyncCounter += (isInSync) ? 0 : 1;
            return ret;
        }
        //END CLASS Demodulator

        //START CLASS complexMath
        bool complexMath::lessThan(std::complex<double> left, std::complex<double> right) {
            return sqrt(left.real() * left.real() + left.imag() * left.imag()) < sqrt(right.real() * right.real() + right.imag() * right.real());
        }
        bool complexMath::biggerThan(std::complex<double> left, std::complex<double> right) {
            return lessThan(right, left);
        }
        //END CLASS complexMath
    }

    namespace decoder {
        //START CLASS UWFinder
        void UWFinder::SetTolerance(int tolerance) {
            this->tolerance = tolerance;
        }
        bool UWFinder::IsFrameDetected(bool lowestFrame, int *nUW, int *rUW, bool *isReversedPolarity, bool *isMidStreamReversePolarity, bool *isReversedFirst) {
            bool isFirstSymbolDetermined = false;
            *isReversedFirst = false;
            *nUW = 0;
            *rUW = 0;
            *isReversedPolarity = false;
            *isMidStreamReversePolarity = false;
            int patternPos = 0;
            // the symbolRegister array stores data as it comes out of the demodulator, 1st symbol at the highest index
            // in the symbolRegister
            // there are 160 + 2 columns per row
            int symbolPos = UWFINDER_FRAME_LENGTH - 1;
            int minPos = 0;
            if (!lowestFrame) {
                symbolPos = 2 * UWFINDER_FRAME_LENGTH - 1;
                minPos = UWFINDER_FRAME_LENGTH;
            }
            for (; symbolPos >= minPos; symbolPos -= 162) {
                //compute normal polarity
                uint8_t pp = nrmPolUwPattern[patternPos];
                uint8_t rpp = revPolUwPattern[patternPos];
                uint8_t sra = symbolRegister[symbolPos];
                uint8_t srb = symbolRegister[symbolPos - 1];
                *nUW += pp ^ sra;
                *nUW += pp ^ srb;
                //compute reverse polarity
                *rUW += rpp ^ sra;
                *rUW += rpp ^ srb;
                if ((!isFirstSymbolDetermined) && (*rUW != *nUW)) {
                    isFirstSymbolDetermined = true;
                    *isReversedFirst = *rUW > *nUW;
                }
                patternPos++;
                //detect polarity change midstream
                if (!*isMidStreamReversePolarity){
                    *isMidStreamReversePolarity = *nUW % 2 != 0;
                }
            }
            *isReversedPolarity = *rUW <= tolerance;
            return (*nUW <= tolerance || *rUW <= tolerance);
        }
        std::vector<UWFinder::uwfinder_result> UWFinder::Decode(uint8_t bitsDemodulated[DEMODULATOR_SYMBOLSPERCHUNK], bool isHardDecision) {
            std::vector<uwfinder_result> ret;
            for(int pos = 0; pos < DEMODULATOR_SYMBOLSPERCHUNK; pos++) {
                //shift register to right with one symbol to make space for the next incoming symbol
                std::copy(symbolRegister, &symbolRegister[(UWFINDER_FRAME_LENGTH*2)], &symbolRegister[1]);
                //store
                symbolRegister[0] = bitsDemodulated[pos];
                symbolCount++;
                if (symbolCount >= UWFINDER_FRAME_LENGTH) {
                    int nUW;
                    int rUW;
                    bool isReversedPolarity;
                    bool isMidStreamReversePolarity;
                    bool isReversedFirst_;
                    if(IsFrameDetected(true, &nUW, &rUW, &isReversedPolarity, &isMidStreamReversePolarity, &isReversedFirst_)) {
                        // We check symbolRegister length
                        // If the length is larger than 1.6 the length of a frame, we can assume
                        // there is another at least 0.6 of a frame before the current detected one
                        // We attempt to correct it and redetect it.
                        // If after correction it passes the tolerence we send it further before the
                        // packet we detected.
                        if (symbolCount > UWFINDER_FRAME_LENGTH * 1.6) {
                            // There is no point in checking if the frame has been detected, however we need the
                            // output values
                            int nUW_;
                            int rUW_;
                            bool isReversedPolarity_;
                            bool isMidStreamReversePolarity_;
                            IsFrameDetected(false, &nUW_, &rUW_, &isReversedPolarity_, &isMidStreamReversePolarity_, &isReversedFirst_);
                            // We reverse all symbols either for reverse polarity or for normal
                            // however we need to know which are which
                            int i_ = 0;
                            if (isReversedFirst_) {
                                // The highest index in the uncertain frame contain reversed polarity symbols
                                i_ = 2 * UWFINDER_FRAME_LENGTH - 81 * rUW_;
                            } else {
                                i_ = 2 * UWFINDER_FRAME_LENGTH - 81 * nUW_;
                            }
                            for (; i_ > UWFINDER_FRAME_LENGTH; i_--) {
                                symbolRegister[i_] = symbolRegister[i_] ^ 0;
                            }
                            // now we reassess the result
                            if(IsFrameDetected(false, &nUW_, &rUW_, &isReversedPolarity_, &isMidStreamReversePolarity_, &isReversedFirst_)) {
                                uwfinder_result res;
                                std::copy(symbolRegister, &symbolRegister[UWFINDER_FRAME_LENGTH*2], res.uwFrame);
                                res.length = UWFINDER_FRAME_LENGTH;
                                res.isReversedPolarity = isReversedPolarity_;
                                //reverse reversed (for hard decision)
                                //this will neede rewriting when the demodulator will output soft symbols.
                                if(isReversedPolarity_) {
                                    for(int i = 0; i < UWFINDER_FRAME_LENGTH; i++) {
                                        res.uwFrame[i] = res.uwFrame[i] ^ 0;
                                    }
                                }
                                res.symbolCount = UWFINDER_FRAME_LENGTH;
                                res.isMidStreamReversePolarity = isMidStreamReversePolarity_;
                                res.BER = std::min(nUW_, rUW_);
                                res.isUncertain = true;
                                res.isHardDecision = isHardDecision;
                                ret.push_back(res);
                            }
                        }
                        uwfinder_result res;
                        for(int i = 0; i < UWFINDER_FRAME_LENGTH; i++) {
                            res.uwFrame[i] = symbolRegister[i];
                        }
                        res.length = UWFINDER_FRAME_LENGTH;
                        res.isReversedPolarity = isReversedPolarity;
                        //reverse reversed (for hard decision)
                        //this will neede rewriting when the demodulator will output soft symbols.
                        if(isReversedPolarity) {
                            for(int i = 0; i < UWFINDER_FRAME_LENGTH; i++) {
                                res.uwFrame[i] = res.uwFrame[i] ^ 0;
                            }
                        }
                        res.symbolCount = symbolCount;
                        res.isMidStreamReversePolarity = isMidStreamReversePolarity;
                        res.BER = std::min(nUW, rUW);
                        res.isUncertain = false;
                        res.isHardDecision = isHardDecision;
                        ret.push_back(res);
                        // There is no point in clearing up the symbolRegister
                        // However also there is no point in detecting until we fill in the register again
                        symbolCount = 0;
                    }
                }
            }
            return ret;
        }
        //END CLASS UWFinder

        //START CLASS Depermuter
        Depermuter::Depermuter() {
            //generate depermutting array
            depermuttingArray.resize(64);
            for (int i = 0; i < 64; i++) {
                //j(i)
                depermuttingArray[i] = (i * 23) % 64;
                //take right-hand 6 bits only
                depermuttingArray[i] = depermuttingArray[i] & 0x3F;
                //get index in source UWFrame (assumes a reversed UWFrame)
                depermuttingArray[i] = depermuttingArray[i] * 162;
            }
        }
        Depermuter::depermuter_result Depermuter::depermute(uint8_t uwFrame[UWFINDER_FRAME_LENGTH], bool isHardDecision) {
            //reverse frame
            std::reverse(uwFrame, &uwFrame[UWFINDER_FRAME_LENGTH]);
            //create the destination
            depermuter_result res;
            //depermute
            for (int i = 0; i < 64; i++) {
                int from_start_index = depermuttingArray[i];
                int to_start_index = (i*162);
                int num_elements = 162;
                for(int k = 0; k < num_elements; k++) {
                    res.depermutedFrame[to_start_index + k] = uwFrame[from_start_index + k];
                }
            }
            res.isHardDecision = isHardDecision;
            res.length = DEPERMUTER_FRAME_LENGTH;

            return res;
        }
        //END CLASS Depermuter

        //START CLASS Deinterleaver
        Deinterleaver::deinterleaver_result Deinterleaver::deinterleave(uint8_t depermutedFrame[DEPERMUTER_FRAME_LENGTH], bool isHardDecision) {
            //store depermutedFrame into the deinterleaver matrix leaving out the UW
            int row = -1;
            int column = 0;
            for (int i = 0; i < DEPERMUTER_FRAME_LENGTH; i++) {
                //at over 162 symbols: reset column, increment row and jump over the UW
                if (i % 162 == 0) {
                    column = 0;
                    row++;
                    i += 2;
                }
                deinterleverMatrix[row][column] = depermutedFrame[i];
                column++;
            }
            //create the destination
            Deinterleaver::deinterleaver_result res;
            //read the matrix into the destination, transposed
            int pos = 0;
            row = 0;
            column = 0;
            for (; row < 64;) {
                res.deinterleavedFrame[pos] = deinterleverMatrix[row][column];
                row++;
                if (row % 64 == 0) {
                    row = 0;
                    column++;
                    if (column == 160) {
                        break;
                    }
                }
                pos++;
            }
            res.isHardDecision = isHardDecision;
            res.length = DEINTERLEAVER_FRAME_LENGTH;
            return res;
        }
        //END CLASS Deinterleaver

        //START CLASS ViterbiDecoder
        ViterbiDecoder::ViterbiDecoder() {
            int amp = 100;
            double esn0 = 5.0;
            double bias = 0.0;
            int scale = 4;
            double noise;
            int s;
            int bit;
            double metrics[2][256];
            double p0;
            double p1;
            // Es/N0 as power ratio
            esn0 = pow(10.0, esn0 / 10);
            noise = 0.5 / esn0; // only half the noise for BPSK
            noise = sqrt(noise); // noise/signal Voltage ratio
            // Zero is a special value, since this sample includes all
            // lower samples that were clipped to this value, i.e., it
            // takes the whole lower tail of the curve
            // P(s|1)
            p1 = (0.5 + 0.5 * std::erf((((0 - offset + 0.5) / amp - 1) / noise) / sqrt2));
            // Prob of this value occurring for a 0-bit
            // P(s|0)
            p0 = (0.5 + 0.5 * std::erf((((0 - offset + 0.5) / amp + 1) / noise) / sqrt2));
            metrics[0][0] = (log(2 * p0 / (p1 + p0)) * log2) - bias;
            metrics[1][0] = (log(2 * p1 / (p1 + p0)) * log2) - bias;
            for (s = 1; s < 255; s++) {
                // P(s|1), prob of receiving s given 1 transmitted
                p1 = (0.5 + 0.5 * std::erf((((s - offset + 0.5) / amp - 1) / noise) / sqrt2))
                        - (0.5 + 0.5 * std::erf((((s - offset - 0.5) / amp - 1) / noise) / sqrt2));
                // P(s|0), prob of receiving s given 0 transmitted
                p0 = (0.5 + 0.5 * std::erf((((s - offset + 0.5) / amp + 1) / noise) / sqrt2))
                        - (0.5 + 0.5 * std::erf((((s - offset - 0.5) / amp + 1) / noise) / sqrt2));
                metrics[0][s] = (log(2 * p0 / (p1 + p0)) * log2) - bias;
                metrics[1][s] = (log(2 * p1 / (p1 + p0)) * log2) - bias;
            }
            // 255 is also a special value
            // P(s|1)
            p1 = 1 - (0.5 + 0.5 * std::erf((((255 - offset - 0.5) / amp - 1) / noise) / sqrt2));
            // P(s|0)
            p0 = 1 - (0.5 + 0.5 * std::erf((((255 - offset - 0.5) / amp + 1) / noise) / sqrt2));
            metrics[0][255] = (log(2 * p0 / (p1 + p0)) * log2) - bias;
            metrics[1][255] = (log(2 * p1 / (p1 + p0)) * log2) - bias;
            // The probability of a raw symbol error is the probability
            // that a 1-bit would be received as a sample with value
            // 0-128. This is the offset normal curve integrated from -Inf to 0.
            for (bit = 0; bit < 2; bit++) {
                for (s = 0; s < 256; s++) {
                    /// Scale and round to nearest integer
                    mettab[bit][s] = (int)floor(metrics[bit][s] * scale + 0.5);
                }
            }
        }
        ViterbiDecoder::viterbidecoder_result ViterbiDecoder::decode(uint8_t deinterleavedFrame[DEINTERLEAVER_FRAME_LENGTH], bool isHardDecision) {
            viterbidecoder_result res;
            uint bitcnt = 0;
            int nbits = DEINTERLEAVER_FRAME_LENGTH / 16;
            int mets[4];
            long bestmetric;
            long beststate;
            int i;
            int j = 0;
            // Initialize arrays
            std::vector<State> state;
            state.resize(64);
            std::vector<State> next;
            next.resize(64);
            std::vector<State> m;
            // Initialize starting metrics to prefered 0 state
            state[0].metric = 0;
            for (i = 1; i < 64; i++) {
                state[i].metric = -999999;
            }
            state[0].path = 0;
            int inputCounter = 0;
            std::vector<uint8_t> input;
            input.resize(DEINTERLEAVER_FRAME_LENGTH);
            /// This is a soft bits viterbi decoder.
            /// If the demodulator outputs hard bits, 0 or 1
            /// we simulate the soft bits as: 28 = zero, 228 = 1
            if (isHardDecision) {
                for (int k = 0; k < DEINTERLEAVER_FRAME_LENGTH; k++) {
                    input[k] = deinterleavedFrame[k] == 0 ? (uint8_t)28 : (uint8_t)228;
                }
            }
            for (bitcnt = 0; bitcnt < nbits * 8; bitcnt++) {
                // Read input symbol pair and compute all possible branch metrics
                mets[0] = mettab[0][input[inputCounter]] + mettab[0][input[inputCounter + 1]];
                mets[1] = mettab[0][input[inputCounter]] + mettab[1][input[inputCounter + 1]];
                mets[2] = mettab[1][input[inputCounter]] + mettab[0][input[inputCounter + 1]];
                mets[3] = mettab[1][input[inputCounter]] + mettab[1][input[inputCounter + 1]];
                inputCounter += 2;
                /// Macro calls originally were generated by genbut.c, as a loop below.
                /// The "C++ to C# Converter" will actually create code for the macros
                /// and we'll end up with about 900+ lines of code, but they will be fully
                /// functional (tested).
                /// Another bonus about using the converter is that it will give you the
                /// indexes into the parity lookup table as well.
                for (i = 0; i < 32; i++) {
                    long m0 = 0;
                    long m1 = 0;
                    int sym = partabIdx[i];
                    /// Add-Compare-Select for 0 branch
                    m0 = state[i].metric + mets[sym];
                    m1 = state[i + 32].metric + mets[3 ^ sym];
                    if (m0 > m1) {
                        next[2 * i].metric = m0;
                        next[2 * i].path = state[i].path << 1;
                    } else {
                        next[2 * i].metric = m1;
                        next[2 * i].path = (state[i + 32].path << 1) | 1;
                    }
                    /// Add-Compare-Select for 1 branch
                    m0 = state[i].metric + mets[3 ^ sym];
                    m1 = state[i + 32].metric + mets[sym];
                    if (m0 > m1) {
                        next[2 * i + 1].metric = m0;
                        next[2 * i + 1].path = state[i].path << 1;
                    } else {
                        next[2 * i + 1].metric = m1;
                        next[2 * i + 1].path = (state[i + 32].path << 1) | 1;
                    }
                }
                // Swap current and next states
                m = state;
                state = next;
                next = m;
                if (bitcnt > DEINTERLEAVER_FRAME_LENGTH - 7) {
                    // In tail, poison non-zero nodes
                    for (i = 1; i < 64; i += 2) {
                        state[i].metric = -9999999;
                    }
                }
                // Produce output every 8 bits once path memory is full
                if ((bitcnt % 8) == 5 && bitcnt > 32) {
                    // Find current best path
                    bestmetric = state[0].metric;
                    beststate = 0;
                    for (i = 1; i < 64; i++) {
                        if (state[i].metric > bestmetric) {
                            bestmetric = state[i].metric;
                            beststate = i;
                        }
                    }
                    res.viterbiFrame[j++] = (uint8_t)(state[beststate].path >> 24);
                }
            }
            // Output remaining bits from 0 state
            if ((i = (int)(bitcnt % 8)) != 6) {
                state[0].path <<= 6 - i;
            }
            res.viterbiFrame[j++] = (uint8_t)(state[0].path >> 24);
            res.viterbiFrame[j++] = (uint8_t)(state[0].path >> 16);
            res.viterbiFrame[j++] = (uint8_t)(state[0].path >> 8);
            res.viterbiFrame[j] = (uint8_t)(state[0].path);
            res.length = VITERBIDECODER_FRAME_LENGTH;
            return res;
        }
        //END CLASS ViterbiDecoder

        //START CLASS Descrambler
        Descrambler::Descrambler() {
            uint8_t x7;
            uint8_t x5;
            uint8_t x4;
            uint8_t x3;
            uint8_t newByte;
            /// Initial state; the documentation found is incorrect as it indicates 0x40, however by
            /// corroboration with more documents and experimenting, the 0x80 generates the correct
            /// scrambling/descrambling array.
            uint8_t register_var = 0x80;
            for (int i = 0; i < DESCRAMBLER_GROUP_COUNT; i++) {
                x7 = (uint8_t)(register_var & 0x01);
                descramblerArray[i] = x7;
                x5 = (uint8_t)((register_var & 0x04) >> 2);
                x4 = (uint8_t)((register_var & 0x08) >> 3);
                x3 = (uint8_t)((register_var & 0x10) >> 4);
                newByte = (uint8_t)(x7 ^ x5 ^ x4 ^ x3);
                register_var >>= 1;
                register_var = (uint8_t)(register_var | (uint8_t)(newByte << 7));
            }
        }
        Descrambler::descrambler_result Descrambler::decode(uint8_t viterbiFrame[VITERBIDECODER_FRAME_LENGTH]) {
            descrambler_result res;


            /// Invert all
            for (int i = 0; i < DESCRAMBLER_FRAME_LENGTH; i++) {
                res.descramblerFrame[i] = invertBits(viterbiFrame[i]);
            }

            /// Apply the bitwise complement only for the "1" states of the descrambler array.
            int j = 0;
            for (int i = 0; i < DESCRAMBLER_GROUP_COUNT; i++) {
                if (descramblerArray[i] == 1) {
                    res.descramblerFrame[j] = (uint8_t)~(res.descramblerFrame[j]);
                    res.descramblerFrame[j + 1] = (uint8_t)~(res.descramblerFrame[j + 1]);
                    res.descramblerFrame[j + 2] = (uint8_t)~(res.descramblerFrame[j + 2]);
                    res.descramblerFrame[j + 3] = (uint8_t)~(res.descramblerFrame[j + 3]);
                }
                j += 4;
            }

            return res;
        }
        uint8_t Descrambler::invertBits(uint8_t input) {
            uint8_t result;
            result = (uint8_t)((input & 0x01));
            result <<= 1;
            result |= (uint8_t)((input & 0x02) >> 1);
            result <<= 1;
            result |= (uint8_t)((input & 0x04) >> 2);
            result <<= 1;
            result |= (uint8_t)((input & 0x08) >> 3);
            result <<= 1;
            result |= (uint8_t)((input & 0x16) >> 4);
            result <<= 1;
            result |= (uint8_t)((input & 0x32) >> 5);
            result <<= 1;
            result |= (uint8_t)((input & 0x64) >> 6);
            result <<= 1;
            result |= (uint8_t)((input & 0x80) >> 7);
            return result;
        }
        //END CLASS Descrambler

        //START CLASS Decoder
        Decoder::Decoder(int tolerance) {
            uwFinder = new UWFinder();
            uwFinder->SetTolerance(tolerance);
            depermuter = new Depermuter();
            deinterleaver = new Deinterleaver();
            viterbiDecoder = new ViterbiDecoder();
            descrambler = new Descrambler();
        }
        std::vector<Descrambler::descrambler_result> Decoder::decode(uint8_t inputBits[DEMODULATOR_SYMBOLSPERCHUNK]) {
            std::vector<Descrambler::descrambler_result> res;
            std::vector<UWFinder::uwfinder_result> uwfinder_result =  uwFinder->Decode(inputBits, false);
            for(int i = 0; i < uwfinder_result.size(); i++) {
                Depermuter::depermuter_result depermuter_result = depermuter->depermute(uwfinder_result[i].uwFrame, true);
                Deinterleaver::deinterleaver_result deinterleaver_result = deinterleaver->deinterleave(depermuter_result.depermutedFrame, true);
                ViterbiDecoder::viterbidecoder_result viterbidecoder_result = viterbiDecoder->decode(deinterleaver_result.deinterleavedFrame, true);
                Descrambler::descrambler_result descrambler_result = descrambler->decode(viterbidecoder_result.viterbiFrame);
                res.push_back(descrambler_result);
            }
            return res;
        }
        //END CLASS Decoder
    }
}
