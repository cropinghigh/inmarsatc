#ifndef INMARSATC_DEMODULATOR_H
#define INMARSATC_DEMODULATOR_H

#if defined(_MSC_VER) || defined(WIN64) || defined(_WIN64) || defined(__WIN64__) || defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#  define Q_DECL_EXPORT __declspec(dllexport)
#  define Q_DECL_IMPORT __declspec(dllimport)
#else
#  define Q_DECL_EXPORT     __attribute__((visibility("default")))
#  define Q_DECL_IMPORT     __attribute__((visibility("default")))
#endif

#if defined(INMARSATC_LIBRARY)
#  define INMARSATC_EXPORT Q_DECL_EXPORT
#else
#  define INMARSATC_EXPORT Q_DECL_IMPORT
#endif

#include <complex>
#include <vector>
#include <algorithm>
#include <math.h>
#include <mutex>

namespace inmarsatc {
    namespace demodulator {
        class FIR {
            public:
                FIR(double b[], int len);
                double filter(double sample);
                double filter();
                void resetFilter();
            private:
                int psLength;
                double* b;
                double* prevSamples;
        };

        class RRC {
            public:
                RRC(double alpha, int firstSize, double sampleRate, double symbolRate);
                double filter(double sample);
            private:
                double* points;
                FIR* rrcFIR;
        };

        class Gardner {
            public:
                Gardner(double sampleRate, double symbolRate);
                bool step(double re, double im, std::complex<double>* output, bool* onPoint);
            private:
                double ts;
                double symbol2xTimer;
                int samplesAgo;
                std::complex<double> maxErrorxAggression;
                int isOnPoint;
                std::complex<double> currentSample;
                std::complex<double> sampleThisOn;
                std::complex<double> error;
                std::complex<double> sampleLastOn;
                std::complex<double> sampleOff;
                double aggression;
                bool result;
                std::complex<double> dummyOutput;
        };

        class CMA {
            public:
                #define CMA_SZ 9
                CMA();
                bool step(std::complex<double> sample, bool isOnPoint, std::complex<double>* output);
                void CMAReset();
                bool result;
            private:
                std::complex<double> dummyOutput;
                std::complex<double> cmaW[CMA_SZ];
                std::complex<double> cmaX[CMA_SZ];
                double beta;
                double stepSize;
                double mean;
                std::complex<double> cmaEqualizerOut;
                std::complex<double> error;
        };

        class AGC {
            public:
                AGC();
                void apply(double* i, double* q);
                void apply(std::complex<double>* value);
                void reset();
            private:
                double mAGC;
        };

        class INMARSATC_EXPORT Demodulator {
            public:
                #define DEMODULATOR_SYMBOLSPERCHUNK 5000
                struct demodulator_result {
                    double meanMagnitude;
                    uint8_t bitsDemodulated[DEMODULATOR_SYMBOLSPERCHUNK];
                };

                Demodulator();
                bool isCmaEnabled();
                bool isAgcEnabled();
                int getLowFreq();
                int getHighFreq();
                double getCenterFreq();
                bool getIsInSync();
                int getNoSyncCount();
                std::complex<double> getScatterPoint();
                void setCmaEnabled(bool cmaEnabled);
                void setAgcEnabled(bool agcEnabled);
                void setLowFreq(int lowFreq);
                void setHighFreq(int highFreq);
                void setCenterFreq(double centerFreq);
                void cmaReset();
                std::vector<demodulator_result> demodulate(std::complex<double> samples[], int length);

            private:
                #define DEMODULATOR_SYMBOLRATE 1200.0
                #define DEMODULATOR_SAMPLERATE 48000.0
                // As the carrier gets locked, the positive and the negative error are equally distributed and
                // their sum gets closer and closer to zero. This is the absolute number I came up with.
                #define DEMODULATOR_MAXABSERRSUM 1.0
                #define DEMODULATOR_ALPHA 0.0065
                // this implementation of a Costas loop tracks easier to the left, so we set the
                // center frequency on purpuse a bit higher that where we expect to find the sync
                #define DEMODULATOR_DEFAULT_FREQ 2600.0
                #define DEMODULATOR_LOW_FREQ 500.0
                #define DEMODULATOR_HIGH_FREQ 4500.0
                double freq = DEMODULATOR_DEFAULT_FREQ;
                std::mutex freq_mtx;
                const double alpha = DEMODULATOR_ALPHA;
                // see literature
                const double beta = alpha * alpha / 4;
                double omega;
                double phase = 0.0;
                int noSyncCounter = 0;
                double error;
                bool isInSync = false;
                int flagcounter;
                std::complex<double> scatterPoint;
                std::complex<double> lastScatterPoint;
                // LPF1, LPF2
                // The literature suggestes that the square-shaped LPF filter should have a
                // bandwidth of no less than half the symbol rate. This will remove the most
                // noise possible without reducing the amplitude of the desired signal at
                // its sampling instants.
                double blpf12[2] = {0.9, 0.9};
                FIR* lpf1;
                FIR* lpf2;
                // The literature indicates that the loop filter should have a
                // response far outside the LPF1/LPF2.The fastest achievable settle
                // time is one in which the VCO has a gain 8x that of the LPF1/LPF2 pole
                // frequency.For LPF3, a factor of six times K or 12 times the LPF1/LPF2
                // pole is a better choice
                double blpf3[2] = {0.48, 0.48};
                FIR* lpf3;
                RRC* rrcRe;
                RRC* rrcIm;
                RRC* rrcRe2;
                RRC* rrcIm2;
                Gardner* gardner;
                CMA* cma;
                bool cmaEnabled = false;
                std::complex<double> gardnerOutputSample;
                std::complex<double> cmaOutputSample;
                uint8_t symbolBuffer[DEMODULATOR_SYMBOLSPERCHUNK];
                AGC* agc;
                bool agcEnabled = false;
                int loFreq = DEMODULATOR_LOW_FREQ;
                int hiFreq = DEMODULATOR_HIGH_FREQ;
        };

        class complexMath {
            public:
                static bool lessThan(std::complex<double> left, std::complex<double> right);
                static bool biggerThan(std::complex<double> left, std::complex<double> right);
        };
    }
}

#endif // INMARSATC_DEMODULATOR_H
