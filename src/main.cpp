#include <iostream>

#include "itkImage.h"
#include "itkImageFileWriter.h"

#include "dce/pksolver.hpp"
#include "core/image.hpp"

using namespace Core::Image;

int main(void) {
	std::cout << "PK Solver!" << std::endl;

	DCE::PKSolver* solver = new DCE::PKSolver(
		"../test/Gd.nii.gz", 
		"../test/preGd.nii.gz", 
		"../test/T1Map.nii.gz", 
		"../test/AIFMask.nii.gz"
	);

	solver->set_time({ 0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 91, 98, 105, 112, 119, 126, 133, 140, 147, 154, 161, 168, 175, 182, 189, 196, 203 });
	solver->set_TR(5.27);
	solver->set_FA(15.0);
	solver->set_relaxivity(0.0037);
	solver->set_model(DCE::TOFTS);

	solver->execute();

	std::cout << "PK Solver: Finished!" << std::endl;

	Image3DPointer KTrans = solver->get_KTrans();
	Image3DPointer Ve = solver->get_Ve();


	using WriterType = itk::ImageFileWriter<Image3DType>;
	itk::SmartPointer<WriterType> writer = WriterType::New();

	writer->SetFileName("../test/output/KTrans.nii.gz");
	writer->SetInput(KTrans);

	try
	{
		writer->Update();
	}
	catch (const itk::ExceptionObject& exp)
	{
		exp.Print(std::cerr);
		return EXIT_FAILURE;  
	}

	itk::SmartPointer<WriterType> writer2 = WriterType::New();

	writer2->SetFileName("../test/output/Ve.nii.gz");
	writer2->SetInput(Ve);

	try
	{
		writer2->Update();
	}
	catch (const itk::ExceptionObject& exp)
	{
		exp.Print(std::cerr);
		return EXIT_FAILURE;  
	}
	return EXIT_SUCCESS; 
}