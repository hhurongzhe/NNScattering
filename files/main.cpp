#include "include/evaluate.hpp"
#include "include/lib_define.hpp"
#include "include/observable.hpp"

int main()
{
    // config file initializing.
    auto ini = inifile_system::inifile("files/NN.ini");
    if (!ini.good())
    {
        std::cerr << ini.error() << std::endl;
        exit(-1);
    }
    auto configs = NN::NN_configs(ini);
    std::cout << configs.result_file() << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    // set parallel threads in openpm.
    constexpr int64_t NUM_THREADS = 16;
    omp_set_num_threads(NUM_THREADS);

    //* main program:

    evaluate::evaluate_phase_np(configs);

    // observable::differential_cross_section(configs);

    // observable::total_cross_section(configs);

    // observable::spin_observables(configs);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout.precision(4);
    std::cout << "Duration: " << duration.count() << " milliseconds\n";
}