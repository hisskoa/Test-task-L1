//	Вариант 1
//1.	Написать на языке С++ класс выполняющий функциональность модулятора QAM(QPSK, QAM16, QAM64)
//2.	Написать на языке С++ класс выполняющий функциональность добавления гауссовского шума к созвездию QAM
//3.	Написать на языке С++ класс выполняющий функциональность демодулятора QAM(QPSK, QAM16, QAM64)
//4.	Написать последовательный вызов 1 - 3 для случайной последовательности бит для разных значений дисперсии шума
//5.	Построить график зависимости вероятности ошибки на бит от  дисперсии шума

#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <fstream>
#include "gnuplot-iostream.h"

class ModulationQPSK {  // Класс выполняющий функциональность модулятора QPSK
private:
    std::vector<double> m_bits;

public:
    ModulationQPSK(std::vector<double> bits) {
        m_bits = bits;
    }
    std::vector<std::complex<double>> Modulation() {
        std::vector<std::complex<double>> qpskBits;

        for (int i{ 0 }; i < m_bits.size() / 2; i++) {
            std::complex<double> z((1 - 2 * m_bits[i * 2]) / sqrt(2), (1 - 2 * m_bits[i * 2 + 1]) / sqrt(2)); // QPSK-модулированный сигнал
            qpskBits.push_back(z);
        }
        return qpskBits;
    }
};
class AddingGaussianNoise {  // Класс выполняющий функциональность добавления гауссовского шума 
public:
    std::vector<std::complex<double>> gaussianNoise(std::vector<std::complex<double>> qpskBits, double snr_dB) {
        double noise_power = pow(10, -0.1 * snr_dB);  // Мощность шума в линейном масштабе 
        std::vector<std::complex<double>> noiseBits;
        const double mean = 0.0;
        const double stddev = noise_power;
        std::default_random_engine generator;
        std::normal_distribution<double> dist(mean, stddev); // Генерирование гауссовского шума

        for (int i{ 0 }; i < qpskBits.size(); i++) {
            std::complex<double> w(dist(generator), dist(generator));
            noiseBits.push_back(qpskBits[i] + w);
        }
        return noiseBits;
    }

};
class DemodulatorQPSK {   // Класс выполняющий функциональность демодулятора
public:
    std::vector<double> Demodulators(std::vector<std::complex<double>> qpskBits) {
        std::vector<double> demodulationSig;
        for (int i{ 0 }; i < qpskBits.size(); i++) {
            if (pow(fabs(qpskBits[i].real() + sqrt(0.5)), 2) - pow(fabs(qpskBits[i].real() - sqrt(0.5)), 2) > 0) {
                demodulationSig.push_back(0);
            }
            else {
                demodulationSig.push_back(1);
            }
            if (pow(fabs(qpskBits[i].imag() + sqrt(0.5)), 2) - pow(fabs(qpskBits[i].imag() - sqrt(0.5)), 2) > 0) {
                demodulationSig.push_back(0);
            }
            else {
                demodulationSig.push_back(1);
            }
        }
        return demodulationSig;
    }

};
int main()
{
    std::random_device rd;
    std::mt19937 gen(rd());  // Затравка для генератора случайных чисел
    std::uniform_int_distribution<int> dist(0, 1);
    int n = 100000; // Количество бит в случайной последовательности
    std::vector<double> input_bits;

    for (int i{ 0 }; i < n; i++) {
        input_bits.push_back(dist(gen));  // Последовательность случайно сгенерированных бит
    }

    ModulationQPSK sig_qpsk(input_bits);
    std::vector<std::complex<double>> modulationBits = sig_qpsk.Modulation();
    int start_snr_dB = 0, end_snr_dB = 12;
    std::vector<int> range_input_snr_dB;
    for (int i{ start_snr_dB }; i < end_snr_dB; i++) {
        range_input_snr_dB.push_back(i);   // Отношение сигнал - шум на поднесущую в децибелах(дБ)
    }
    std::vector<double> bitErrorRate;

    for (int i{ 0 }; i < range_input_snr_dB.size(); i++) {
        double input_snr_dB = range_input_snr_dB[i];
        double count_bit_error = 0;
        AddingGaussianNoise sig_noise;
        std::vector<std::complex<double>> noiseBits = sig_noise.gaussianNoise(modulationBits, input_snr_dB);

        DemodulatorQPSK llr;
        std::vector<double> demodulationBits = llr.Demodulators(noiseBits);

        for (int j{ 0 }; j < input_bits.size(); j++) {
            if (input_bits[j] != demodulationBits[j]) {
                count_bit_error = count_bit_error + 1;
            }
        }
        bitErrorRate.push_back(count_bit_error / input_bits.size());
    }


    // Сохраняем данные в файл
    std::ofstream fout("Output data.txt");
    if (fout.is_open()) {
        for (int i{ 0 }; i < bitErrorRate.size(); i++)
            fout << range_input_snr_dB[i] << "\t" << bitErrorRate[i] << std::endl;
    }
    fout.close();

    // Строим график
    Gnuplot gp("\"bin\\gnuplot.exe\"");

    gp << "set logscale y\n";
    gp << "set ylabel 'BER'\n";
    gp << "set xlabel 'SNR'\n";
    gp << "plot 'Output data.txt' u 1:2 w linesp" << std::endl;
    std::cin.get();

    return 0;
}
