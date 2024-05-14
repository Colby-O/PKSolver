#include "core/image.hpp"
    
namespace Core::Image {
    void image_4d_append_to_front(const Image3DPointer& image3D, const Image4DPointer& image4D, Image4DPointer& resultImage) {
        // Get the size of the original 4D image
        Image4DType::SizeType size4D = image4D->GetLargestPossibleRegion().GetSize();

        // Create a new 4D image with size adjusted to accommodate the 3D image
        Image4DType::SizeType newSize;
        newSize[0] = size4D[0];
        newSize[1] = size4D[1];
        newSize[2] = size4D[2];
        newSize[3] = size4D[3] + 1; // Increase the size of the fourth dimension
        Image4DType::RegionType newRegion;
        newRegion.SetSize(newSize);
        resultImage->SetRegions(newRegion);
        resultImage->Allocate();

        // Copy data from the 3D image to the start of the 4D image
        itk::ImageRegionConstIterator<Image3DType> inputIt(image3D, image3D->GetLargestPossibleRegion());
        itk::ImageRegionIterator<Image4DType> outputIt(resultImage, resultImage->GetLargestPossibleRegion());

        while (!inputIt.IsAtEnd()) {
            outputIt.Set(inputIt.Get());
            ++inputIt;
            ++outputIt;
        }

        // Copy data from the original 4D image to the new 4D image starting from the second position in the fourth dimension
        Image4DType::IndexType index4D;
        index4D.Fill(0);
        index4D[3] = 1; // Start from the second position in the fourth dimension

        itk::ImageRegionIterator<Image4DType> originalIt(image4D, image4D->GetLargestPossibleRegion());
        itk::ImageRegionIterator<Image4DType> resultIt(resultImage, resultImage->GetLargestPossibleRegion());
        resultIt.SetIndex(index4D); // Set the iterator to start from the second position in the fourth dimension

        while (!originalIt.IsAtEnd()) {
            resultIt.Set(originalIt.Get());
            ++originalIt;
            ++resultIt;
        }
    }

    void image_4d_mean(const Image4DPointer& image4D, Image3DPointer& resultImage) {
        // Get the size of the original 4D image
        Image4DType::SizeType size4D = image4D->GetLargestPossibleRegion().GetSize();

        // Create a new 3D image with the same size as the first three dimensions of the original 4D image
        Image3DType::SizeType size3D;
        size3D[0] = size4D[0];
        size3D[1] = size4D[1];
        size3D[2] = size4D[2];
        Image3DType::RegionType region;
        region.SetSize(size3D);
        resultImage->SetRegions(region);
        resultImage->Allocate();

        // Calculate the mean along the last dimension for each voxel in the resulting 3D image
        itk::ImageRegionIteratorWithIndex<Image3DType> outputIt(resultImage, resultImage->GetLargestPossibleRegion());
        while (!outputIt.IsAtEnd()) {
            // Get the index of the current voxel
            Image3DType::IndexType index3D = outputIt.GetIndex();
            Image4DType::IndexType index = Image4DType::IndexType();
            index[0] = index3D[0];
            index[1] = index3D[1];
            index[2] = index3D[2];

            // Calculate the mean value along the last dimension for the current voxel
            double sum = 0.0;
            for (unsigned int i = 0; i < size4D[3]; ++i) {
                index[3] = i;
                sum += image4D->GetPixel(index);
            }
            double mean = sum / size4D[3];

            // Set the mean value to the current voxel in the resulting 3D image
            outputIt.Set(mean);

            ++outputIt;
        }
    }
}