#include "dce/pksolver.hpp"

#define IS_NAN(x) ((x) != (x))
#define M_PI 3.14159265358979323846

using namespace Core::Image;

std::vector<double> DCE::sample_function(double(*function)(double, double, double), double startTime, double endTime, double samplingRate, double Ktrans, double Ve) {
    std::vector<double> samples;
    int numSamples = std::ceil((endTime - startTime) / samplingRate);
    samples.reserve(numSamples);
    for (int i = 0; i <= numSamples; ++i) {
        double t = startTime + i * samplingRate;
        samples.push_back(function(t, Ktrans, Ve));
    }
    return samples;
}

double DCE::tissue_impulse_response(double t, double Ktrans, double Ve) {
    return Ktrans * std::exp(-(Ktrans / Ve) * t);
}

std::vector<double> DCE::fft_convolution(std::vector<double> AIF, std::vector<double> tissueImpulseResponse) {
    int N = AIF.size();

    // Allocate memory for FFT inputs and outputs
    fftw_complex* AIF_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* tissueImpulseResponse_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex* result_fft = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * N);

    // Create FFTW plans
    fftw_plan AIF_plan = fftw_plan_dft_r2c_1d(N, AIF.data(), AIF_fft, FFTW_ESTIMATE);
    fftw_plan tissueImpulseResponse_plan = fftw_plan_dft_r2c_1d(N, tissueImpulseResponse.data(), tissueImpulseResponse_fft, FFTW_ESTIMATE);

    // Execute FFT plans
    fftw_execute(AIF_plan);
    fftw_execute(tissueImpulseResponse_plan);

    // Element-wise multiplication in frequency domain
    for (size_t i = 0; i < N; ++i) {
        result_fft[i][0] = AIF_fft[i][0] * tissueImpulseResponse_fft[i][0] - AIF_fft[i][1] * tissueImpulseResponse_fft[i][1];
        result_fft[i][1] = AIF_fft[i][0] * tissueImpulseResponse_fft[i][1] + AIF_fft[i][1] * tissueImpulseResponse_fft[i][0];
    }

    // Allocate memory for the result of the IFFT
    double* result = (double*)fftw_malloc(sizeof(double) * N);

    // Create IFFT plan
    fftw_plan result_plan = fftw_plan_dft_c2r_1d(N, result_fft, result, FFTW_ESTIMATE);

    // Execute IFFT plan
    fftw_execute(result_plan);

    // Copy the real part of the result
    std::vector<double> convolutionResult(N);
    for (size_t i = 0; i < N; ++i) {
        convolutionResult[i] = result[i] / N;
    }

    // Free memory and destroy plans
    fftw_destroy_plan(AIF_plan);
    fftw_destroy_plan(tissueImpulseResponse_plan);
    fftw_destroy_plan(result_plan);
    fftw_free(AIF_fft);
    fftw_free(tissueImpulseResponse_fft);
    fftw_free(result_fft);
    fftw_free(result);

    return convolutionResult;
}

DCE::PKParameters::PKParameters() {

}

DCE::PKParameters::PKParameters(DCE::PKModelType modelType) : DCE::PKParameters::PKParameters() {
    set_model(modelType);
}

void DCE::PKParameters::set_model(DCE::PKModelType modelType) {
    m_ModelType = modelType;
}

DCE::PKModelType DCE::PKParameters::get_model() {
    return m_ModelType;
}

void DCE::PKParameters::set_parameters(const std::vector<double>& params) {
    switch (m_ModelType) {
        case TOFTS:
            T_Ktrans = params[0];
            T_Kep = params[1];
            T_ve = params[2];
            break;
        case EXTENDED_TOFTS:
            ET_Ktrans = params[0];
            ET_Kep = params[1];
            ET_ve = params[2];
            ET_vp = params[3];
            break;
        default:
            break;
    }
}

std::vector<double> DCE::PKParameters::get_parameters() {
    switch (m_ModelType) {
        case TOFTS:
            return std::vector<double>({T_Ktrans, T_Kep, T_ve});
        case EXTENDED_TOFTS:
            return std::vector<double>({ET_Ktrans, ET_Kep, ET_ve, ET_vp});
        default:
            return std::vector<double>();
    }
}


DCE::PKSolver::PKSolver(
    const std::string& GdFilename, 
    const std::string& preGdFilename, 
    const std::string& T10Filename, 
    const std::string& AIFMaskFilename
) {
    // Loads pre-contrast  T1 values
    m_T10 = Core::Image::load<double, 3>(T10Filename);
    // Loads AIF Mask
    m_AIFMask  = Core::Image::load<std::uint8_t, 3>(AIFMaskFilename);
    // Loads N Signal measurement pre-contrast 
	Image4DPointer S0_raw = Core::Image::load<double, 4>(preGdFilename);
    // Loads Signal measurement post-contrast 
	Image4DPointer signal  = Core::Image::load<double, 4>(GdFilename);

    // Calculates the means signal pre-contrast of each of the measurements
    // If there is more than one.
	Image3DPointer S0 = Image3DType::New();
	image_4d_mean(S0_raw , S0);

    // Appends the pre-contrast to the front of the post-contrast singal volume.
	m_S = Image4DType::New();
	image_4d_append_to_front(S0, signal, m_S);
}

void DCE::PKSolver::set_TR(double TR) {
    m_TR = TR;
}

void DCE::PKSolver::set_FA(double FA){
    m_FA = FA;
}

void DCE::PKSolver::set_relaxivity(double relaxivity){
    m_relaxivity = relaxivity;
}

void DCE::PKSolver::set_model(PKModelType modelType) {
    m_ModelType = modelType;
}

void DCE::PKSolver::set_time(std::vector<double> time) {
    m_Time = time;
}

Image3DPointer DCE::PKSolver::get_KTrans() {
    return m_KTrans;
}

Image3DPointer DCE::PKSolver::get_Ve() {
    return m_Ve;
}

Image3DPointer DCE::PKSolver::get_Kep() {
    return m_Kep;
}

void DCE::PKSolver::execute() {
    // Initalizing concentration map.
    Image4DType::SizeType size = m_S->GetLargestPossibleRegion().GetSize();

    calculate_concentration_time_curves(m_TR, m_FA, m_relaxivity);
    
    // Calcautes the AIF from an averge of the concentration time curves in the AIFMask.
    std::vector<double> Cp = calculate_AIF();

    m_KTrans = create_image<double, 3>(m_T10->GetLargestPossibleRegion().GetSize());
    m_Kep = create_image<double, 3>(m_T10->GetLargestPossibleRegion().GetSize());
    m_Ve = create_image<double, 3>(m_T10->GetLargestPossibleRegion().GetSize());

    //for (int i = 0; i < Cp.size(); i++) std::cout << "AIF at timepoint " << i << " is " << Cp[i] << std::endl; 

    for (unsigned int i = 0; i < size[0]; ++i) {
        std::cout << i << "/" << size[0] << std::endl;
        for (unsigned int j = 0; j < size[1]; ++j) {
            for (unsigned int k = 0; k < size[2]; ++k) {
                std::vector<double> Cb = std::vector<double>();
                for (unsigned int l = 0; l < size[3] - 1; ++l) {
                    Cb.push_back(m_ConcentrationMap->GetPixel({i, j, k, l}));
                }

                DCE::PKParameters* params = fit(Cb, Cp, m_Time);

                //std::cout << "KTrans " << params->get_parameters()[0] << " Ve " << params->get_parameters()[1] << std::endl; 

                m_KTrans->SetPixel({i, j, k}, params->get_parameters()[0]);
                m_Ve->SetPixel({i, j, k}, params->get_parameters()[1]);
                m_Kep->SetPixel({i, j, k}, params->get_parameters()[2]);
            }
        }
    }
}

DCE::PKParameters* DCE::PKSolver::fit(
    const std::vector<double>& Cb,
    const std::vector<double>& Cp,
    const std::vector<double>& time
) {
    assert(m_Time.size() == Cb.size() && m_Time.size() == Cp.size());

    // Define the cost function for optimization
    CostFunction* costFunction = new CostFunction(Cb);
    costFunction->set_AIF(Cp);
    costFunction->set_sampling_parameters(m_Time[0], m_Time[time.size() - 1], m_Time[1] - m_Time[0]);

    // Set up Levenberg-Marquardt optimizer
    itk::LevenbergMarquardtOptimizer::Pointer optimizer = itk::LevenbergMarquardtOptimizer::New();

    // Set initial parameters and parameter scales
    itk::LevenbergMarquardtOptimizer::ParametersType initialParameters(2);

    // Inital Guess KTrans
    initialParameters[0] = 0.1;
    // Inital Guess Ve
    initialParameters[1] = 0.2;

    //itk::OptimizerParameterScalesEstimator::Pointer scalesEstimator = itk::OptimizerParameterScalesEstimator::New();
    //scalesEstimator->SetMetricSamplingPercentage(1.0); // Use all the data points for scale estimation
    //scalesEstimator->SetMetricSamplingPercentage(10); // Use 10% of the data points for estimating the scales
    //optimizer->SetScalesEstimator(scalesEstimator);

    optimizer->SetValueTolerance(1e-4); 
    optimizer->SetGradientTolerance(1e-4);
    //optimizer->set_g_tolerance(1e-4); 
    //optimizer->set_x_tolerance(1e-5); 
    optimizer->SetEpsilonFunction(1e-9); 
    optimizer->SetNumberOfIterations(100);

    // Set other optimizer parameters
    optimizer->SetInitialPosition(initialParameters);
    optimizer->UseCostFunctionGradientOff();
    optimizer->SetCostFunction(costFunction);
    //optimizer->SetNumberOfIterations(100);

    // Perform optimization
    try {
        optimizer->StartOptimization();
    } catch (itk::ExceptionObject& error) {
        std::cerr << "Error: " << error << std::endl;
        return nullptr;
    }

    itk::LevenbergMarquardtOptimizer::ParametersType optimizedParameters = optimizer->GetCurrentPosition();
    PKParameters* pkParameters = new PKParameters();
    pkParameters->set_model(m_ModelType);
    pkParameters->set_parameters({optimizedParameters[0], optimizedParameters[1], optimizedParameters[0] / optimizedParameters[1]});
    return pkParameters;
}

void DCE::PKSolver::calculate_concentration_time_curves(
    double TR,
    double flip,
    double relaxivity
) {
    Image4DType::SizeType size = m_S->GetLargestPossibleRegion().GetSize();;
    
    size[3]--;
    m_ConcentrationMap = create_image<double, 4>(size);
    size[3]++;
    
    for (unsigned int i = 0; i < size[0]; ++i) {
        for (unsigned int j = 0; j < size[1]; ++j) {
            for (unsigned int k = 0; k < size[2]; ++k) {
                for (unsigned int l = 1; l < size[3]; ++l) {
                    // Calculating the concentration of Gd in blood.
                    double C = convert_signal_to_concentration(
                        m_S->GetPixel({{i, j, k, l}}), 
                        m_S->GetPixel({{i, j, k, 0}}),
                        m_T10->GetPixel({{i, j, k}}),
                        TR,
                        flip,
                        relaxivity
                    );
                    m_ConcentrationMap->SetPixel({i, j, k, l - 1}, C);
                }
            }
        }
    }
}

std::vector<double> DCE::PKSolver::calculate_AIF() {
    Image4DType::SizeType size = m_S->GetLargestPossibleRegion().GetSize();;

    std::vector<double> Cp = std::vector<double>();
    for (unsigned int i = 0; i < size[3] - 1; ++i) Cp.push_back(0.0);

    int numberPixels = 0;

    for (unsigned int i = 0; i < size[0]; ++i) {
        for (unsigned int j = 0; j < size[1]; ++j) {
            for (unsigned int k = 0; k < size[2]; ++k) {
                if (m_AIFMask->GetPixel({{i, j, k}}) == 0) continue;
                numberPixels++;
                for (unsigned int l = 0; l < size[3] - 1; ++l) {
                    // Calculating the concentration of Gd in blood.
                    double C = m_ConcentrationMap->GetPixel({i, j, k, l});
                    Cp[l] += C;
                }
            }
        }
    }

    for (unsigned int i = 0; i < size[3] - 1; ++i) Cp[i] /= numberPixels;

    return Cp;
}

Image4DPointer DCE::PKSolver::get_concentration_map() {
    return m_ConcentrationMap;
}

std::vector<double> DCE::PKSolver::compute_derivative() {

}

double DCE::PKSolver::calulate_s0() {

}

/*
*    Derives the concentration in blood of Gd using Spoiled Gradient Echo Sequences Equation 
*    Ignoring T2* effects and assuming T2* >>> TE T1(t) can be given by the following
*       
*       E10 = exp(-TR / T10)
*       B = (1 - E10) / (1 - cos(FA) * E10)
*       A = B * S(t) / S(0)
*       1/T1(t) = (-1 / TR) * ln((1 - A) / (1 - cos(FA) * A))
*    
*    Finally the concentration at time t can be computed as follows:
*       
*       C(t) = (1/T1(t) - 1/T(0)) / RGd

*    Where RGd in the relaxivity of Gd (obtained from contrast agent manufacturer specifications).
*    
*    @param S          - Signal at time t
*    @param S0         - Signal before concentration agent injection
*    @param T10        - T1 before concentration agent injection
*    @param TR         - Repetition time in ms
*    @param flip       - Flip angle in degrees
*    @param relaxivity - Constrast agent relaxivity
*    @return concentration in blood of Gd at specific timepoint.
*/
double DCE::PKSolver::convert_signal_to_concentration(
    double S, 
    double S0, 
    double T10,
    double TR,
    double flip,
    double relaxivity
) {
    // Checking for divide by zero
    if (S0 == 0 || T10 == 0) {
        return 0;
    }

    // Converting the flip angle from degrees to radians
    double flipRad = flip * (M_PI / 180.0);

    // Calcauted parameters used to computed T1(t)
    double E10 = std::exp(-TR / T10);
    double cosFA = std::cos(flipRad);
    double B = (1 - E10) / (1 - cosFA * E10);
    double A = B * S / S0;

    double log = std::log((1 - A) / (1 - cosFA * A));

    // Check for nan resulting from the above log
    if (IS_NAN(log)) {
        return 0;
    }

    // Computed the concentration of Gd in blood.
    double recipT10 = 1 / T10;
    double recipT1 = (-1.0 / TR) * log;
    double C = (recipT1 - recipT10) / relaxivity;

    // Ensures the measured concentration is not null
    assert(!IS_NAN(C));

    // Checks if the concentration is neagtive and if so return zero
    return (C < 0) ? 0.0 : C;
}