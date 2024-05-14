#ifndef IMAGE_HPP_
#define IMAGE_HPP_

#include <string>
#include <cstdint>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"

namespace Core::Image {
	using BinaryMask = itk::Image<std::uint8_t, 3>;
	using Image3DType = itk::Image<double, 3>;
	using Image4DType = itk::Image<double, 4>;
	using BinaryMaskPointer = itk::SmartPointer<BinaryMask>;
	using Image3DPointer = itk::SmartPointer<Image3DType>;
	using Image4DPointer = itk::SmartPointer<Image4DType>;

	/*
	*    Reads an image in as an itk::Image
	*    
	*    @param filename to be loaded
	*    @return ITK image of type itk::Image<T, dim>
	*/
	template<typename T, unsigned int dim, typename ImageType = itk::Image<T, dim>>
	itk::SmartPointer<ImageType> load(const std::string& filename) {
		// Create an image reader
		using ReaderType = itk::ImageFileReader<ImageType>;
		auto reader = ReaderType::New();
		// Reads the image
		reader->SetFileName(filename.c_str());
		reader->Update();

		return reader->GetOutput();
	}

	/*
	*    Resmaples an image (src) to a reference image (ref) using a linear interpolator.
	*    
	*    @param srouce image (to be resampled) and the reference image.
	*    @return The resampled source image 
	*/
	template<typename T, unsigned int dim, typename ImageType = itk::Image<T, dim>, typename Pointer = itk::SmartPointer<ImageType>>
	Pointer resample(Pointer src, Pointer ref) {
		// Type defs for different elements of resampling
		using ResampleFilterType = itk::ResampleImageFilter<ImageType, ImageType>;
		using InterpolatorType = itk::LinearInterpolateImageFunction<ImageType, double>;
		using TransformType = itk::IdentityTransform<double, 2>;

		// Preparing the resamplier
		//auto transform = TransformType::New();
		//transform->SetIdentity();
		
		auto interpolator = InterpolatorType::New();
		
		auto resampler = ResampleFilterType::New();
		//resampler->SetTransform(&*transform);
		resampler->SetInterpolator(interpolator);
		resampler->SetInput(src);
		resampler->SetUseReferenceImage(true);
		resampler->SetReferenceImage(ref);
		resampler->SetDefaultPixelValue(0);
		resampler->Update();
		
		// Returns the resampled image
		return resampler->GetOutput();
	}

	template<typename T, unsigned int dim, typename ImageType = itk::Image<T, dim>, typename Pointer = itk::SmartPointer<ImageType>>
	Pointer create_image(typename ImageType::SizeType size) {
		Pointer img = ImageType::New();
		img->SetRegions(size); 
		img->Allocate();
		return img;
	}

	void image_4d_append_to_front(const Image3DPointer& image3D, const Image4DPointer& image4D, Image4DPointer& resultImage);
	void image_4d_mean(const Image4DPointer& image4D, Image3DPointer& resultImage);
}
#endif
