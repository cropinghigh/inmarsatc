#include "inmarsatc.h"

namespace inmarsatc {
    namespace demodulator {
        //START CLASS FIR
        FIR::FIR(double b[], int len) {
            this->b = new double[len];
            std::copy(b, &b[len-1], this->b);
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
                // keep frequency in range(round)
                if(freq < loFreq) {
                    freq = hiFreq;
                    omega = 2.0 * M_PI * freq / DEMODULATOR_SAMPLERATE;
                }
                if(freq > hiFreq) {
                    freq = loFreq;
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
                uint8_t iBit = (((scatterPoint.real() > 0.0) - (scatterPoint.real() < 0.0)) + 1) / 2;
                //uint8_t iBit = (((scatterPoint.real() > 0.0) - (scatterPoint.real() < 0.0)) + 1.0) / 2.0;
                //rotate right
                std::copy(symbolBuffer, &symbolBuffer[DEMODULATOR_SYMBOLSPERCHUNK-1], &symbolBuffer[1]);
                symbolBuffer[0] = iBit;
                flagcounter++;
                if(flagcounter == DEMODULATOR_SYMBOLSPERCHUNK) {
                    demodulator_result res;
                    std::copy(symbolBuffer, &symbolBuffer[DEMODULATOR_SYMBOLSPERCHUNK-1], res.bitsDemodulated);
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
                std::copy(symbolRegister, &symbolRegister[(UWFINDER_FRAME_LENGTH*2)-1], &symbolRegister[1]);
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
                                symbolRegister[i_] = symbolRegister[i_] ^ 1;
                            }
                            // now we reassess the result
                            if(IsFrameDetected(false, &nUW_, &rUW_, &isReversedPolarity_, &isMidStreamReversePolarity_, &isReversedFirst_)) {
                                uwfinder_result res;
                                std::copy(symbolRegister, &symbolRegister[UWFINDER_FRAME_LENGTH-1], res.uwFrame);
                                res.length = UWFINDER_FRAME_LENGTH;
                                res.isReversedPolarity = isReversedPolarity_;
                                //reverse reversed (for hard decision)
                                //this will neede rewriting when the demodulator will output soft symbols.
                                if(isReversedPolarity_) {
                                    for(int i = 0; i < UWFINDER_FRAME_LENGTH; i++) {
                                        res.uwFrame[i] = res.uwFrame[i] ^ 1;
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
                                res.uwFrame[i] = res.uwFrame[i] ^ 1;
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
            int bitcnt = 0;
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
            res.frameNumber = res.descramblerFrame[2] << 8 | res.descramblerFrame[3];
            res.timestamp = std::chrono::system_clock::now();
            res.length = DESCRAMBLER_FRAME_LENGTH;
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
        std::vector<Decoder::decoder_result> Decoder::decode(uint8_t inputBits[DEMODULATOR_SYMBOLSPERCHUNK]) {
            std::vector<decoder_result> ret;
            std::vector<UWFinder::uwfinder_result> uwfinder_result =  uwFinder->Decode(inputBits, false);
            for(int i = 0; i < (int)uwfinder_result.size(); i++) {
                Depermuter::depermuter_result depermuter_result = depermuter->depermute(uwfinder_result[i].uwFrame, true);
                Deinterleaver::deinterleaver_result deinterleaver_result = deinterleaver->deinterleave(depermuter_result.depermutedFrame, true);
                ViterbiDecoder::viterbidecoder_result viterbidecoder_result = viterbiDecoder->decode(deinterleaver_result.deinterleavedFrame, true);
                Descrambler::descrambler_result descrambler_result = descrambler->decode(viterbidecoder_result.viterbiFrame);
                decoder_result res;
                std::copy(descrambler_result.descramblerFrame, &descrambler_result.descramblerFrame[DESCRAMBLER_FRAME_LENGTH-1], res.decodedFrame);
                res.length = descrambler_result.length;
                res.BER = uwfinder_result[i].BER;
                res.frameNumber = descrambler_result.frameNumber;
                res.isHardDecision = true;
                res.isMidStreamReversePolarity = uwfinder_result[i].isMidStreamReversePolarity;
                res.isReversedPolarity = uwfinder_result[i].isReversedPolarity;
                res.isUncertain = uwfinder_result[i].isUncertain;
                res.timestamp = descrambler_result.timestamp;
                ret.push_back(res);
            }
            return ret;
        }
        //END CLASS Decoder
    }

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
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << inputFrame.decodedFrame[*pos + 6 + i];
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
            double uplinkChannelMhz = ((inputFrame.decodedFrame[*pos + 3] << 8 | inputFrame.decodedFrame[*pos + 3 + 1]) - 6000) * 0.0025 + 1626.5;
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
            double uplinkChannelMhz = ((inputFrame.decodedFrame[*pos + 2] << 8 | inputFrame.decodedFrame[*pos + 2 + 1]) - 6000) * 0.0025 + 1626.5;
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
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << inputFrame.decodedFrame[*pos + 8 + i];
            }
            std::string unknown1Hex = os.str();
            os.clear();
            for(int i = 0; i < 4; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << inputFrame.decodedFrame[*pos + 10 + i];
            }
            std::string unknown2Hex = os.str();
            os.clear();
            for(int i = 0; i < 2; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << inputFrame.decodedFrame[*pos + 15 + i];
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
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << inputFrame.decodedFrame[*pos + 2 + i];
            }
            std::string les = os.str();
            double downlinkChannelMhz = ((inputFrame.decodedFrame[*pos + 5] << 8 | inputFrame.decodedFrame[*pos + 5 + 1]) - 8000) * 0.0025 + 1530.5;
            os.clear();
            for(int i = 0; i < 1; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << inputFrame.decodedFrame[*pos + 7 + i];
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
                    os << std::setfill('0') << std::setw(2) << std::right << std::hex << inputFrame.decodedFrame[*pos + 6 + i];
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
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << inputFrame.decodedFrame[*pos + 6 + i];
            }
            std::string unknown1Hex = os.str();
            int shortMessageLength = inputFrame.decodedFrame[*pos + 9];
            os.clear();
            for(int i = 0; i < 1; i++) {
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << inputFrame.decodedFrame[*pos + 10 + i];
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
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << inputFrame.decodedFrame[*pos + 2 + i];
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
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << address[i + 1];
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
                os << std::setfill('0') << std::setw(2) << std::right << std::hex << address[i + 1];
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
                        case 0x03:
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
