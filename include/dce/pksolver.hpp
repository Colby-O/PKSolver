#ifndef PKSOLVER_HPP_
#define PKSOLVER_HPP_

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fftw3.h>

#include "itkImage.h"
#include "itkLevenbergMarquardtOptimizer.h"

#include "core/image.hpp"

namespace DCE 
{
    using namespace Core::Image;

    
    // Sample the continuous-time function at regular intervals
    std::vector<double> sample_function(double(*function)(double, double, double), double startTime, double endTime, double samplingRate, double Ktrans, double Ve);

    // Define the continuous-time tissue impulse response function with Ve
    double tissue_impulse_response(double t, double Ktrans, double Ve);

    // Perform FFT-based convolution
    std::vector<double> fft_convolution(std::vector<double> AIF, std::vector<double> tissueImpulseResponse);
    
    enum PKModelType {
        TOFTS,
        EXTENDED_TOFTS
    };

    class PKParameters {
    private:
        PKModelType m_ModelType;  

        std::vector<double> m_Cb;
        std::vector<double> m_Cp;
        std::vector<double> m_time;

        union {
            struct {
                double ET_Ktrans, ET_Kep, ET_ve, ET_vp;
            };
            struct {
                double T_Ktrans, T_Kep, T_ve;
            };
        };
    public:
        PKParameters();
        PKParameters(PKModelType modelType);

        void set_model(PKModelType modelType);
        PKModelType get_model();
        void set_parameters(const std::vector<double>& params);
        std::vector<double> get_parameters();
    };
    
    class CostFunction : public itk::MultipleValuedCostFunction {
    public:
        using Self = CostFunction;
        using Superclass = itk::MultipleValuedCostFunction;
        using Pointer = itk::SmartPointer<Self>;
        using ParametersType = Superclass::ParametersType;
        using MeasureType = Superclass::MeasureType;
        using DerivativeType = Superclass::DerivativeType;

        CostFunction(const std::vector<double>& observedData) : m_ObservedData(observedData) { }

        MeasureType GetValue(const ParametersType & parameters) const override {

            std::vector<double> predictedData = fft_convolution(m_AIF, sample_function(tissue_impulse_response, m_StartTime, m_EndTime, m_SamplingRate, parameters[0], parameters[1]));

            MeasureType result(GetNumberOfValues());
            for (size_t i = 0; i < m_ObservedData.size(); ++i) {
                //std::cout << "Pred: " << predictedData[i] << " " << parameters[0] << " " << parameters[1] << std::endl;
                //std::cout<< "Cost: "  << std::pow(m_ObservedData[i] - predictedData[i], 2) << std::endl;
                result[i] = std::pow(m_ObservedData[i] - predictedData[i], 2);
            }

            return result;
        }

        void GetDerivative(const ParametersType & parameters, DerivativeType & derivative) const override {

        }

        unsigned int GetNumberOfParameters() const override {
            return 2;
        }

        unsigned int GetNumberOfValues() const override {
            return m_ObservedData.size();
        }

        void set_AIF(const std::vector<double>& AIF) { m_AIF = AIF; }

        void set_sampling_parameters(double startTime, double endTime, double samplingRate) { 
            m_StartTime = startTime; 
            m_EndTime = endTime; 
            m_SamplingRate = samplingRate; 
        }

    private:
        const std::vector<double>& m_ObservedData;
        std::vector<double> m_AIF;
        double m_StartTime;
        double m_EndTime;
        double m_SamplingRate;
    };
    
    class PKSolver 
    {
    public:
        PKSolver(
            const std::string& GdFilename, 
            const std::string& preGdFilename, 
            const std::string& T10Filename, 
            const std::string& AIFMaskFilename
        );

        void execute();
        Image4DPointer get_concentration_map();
        Image3DPointer get_KTrans();
        Image3DPointer get_Ve();
        Image3DPointer get_Kep();
        void set_time(std::vector<double> time);
        void set_TR(double TR);
        void set_FA(double FA);
        void set_relaxivity(double relaxivity);
        void set_model(PKModelType modelType);
    private:
        PKModelType m_ModelType;
        Image4DPointer m_S, m_ConcentrationMap; 
        Image3DPointer m_T10, m_KTrans, m_Ve, m_Kep; 
        BinaryMaskPointer m_AIFMask;

        double m_TR, m_FA, m_relaxivity;
        std::vector<double> m_Time;

        double compute_bolus_arrive_time();

        double calulate_s0();

        void calculate_concentration_time_curves(
            double TR,
            double flip,
            double relaxivity
        );

        double convert_signal_to_concentration(
            double S, 
            double S0, 
            double T1,
            double TR,
            double flip,
            double relaxivity
        );

        PKParameters* fit(
            const std::vector<double>& Cb,
            const std::vector<double>& Cp,
            const std::vector<double>& time
        );

        std::vector<double> calculate_AIF();
    };
}
#endif