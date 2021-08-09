#include <inmarsatc_demodulator.h>

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
            //this->freq = centerFreq;
            this->omega = (centerFreq * (2.0 * M_PI)) / DEMODULATOR_SAMPLERATE;
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
                meanMagnitude = 1; //avoid zero division
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
                phase += omega + alpha * error;
                // carrier
                omega += beta * error;
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
                vI = cos(-phase);
                vQ = sin(-phase);
                // mixer
                I = vI * samples[i].real();
                Q = vQ * samples[i].imag();
                // LPFs
                I = lpf1->filter(I);
                Q = lpf2->filter(Q);
                // VCO error
                //error = I * Q;
                error = tanh(I * Q);
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
}
