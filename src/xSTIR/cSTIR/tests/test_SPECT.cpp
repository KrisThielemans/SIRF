/*
SyneRBI Synergistic Image Reconstruction Framework (SIRF)
Copyright 2018 - 2020 Rutherford Appleton Laboratory STFC
Copyright 2018 - 2020 University College London

This is software developed for the Collaborative Computational
Project in Synergistic Reconstruction for Biomedical Imaging (formerly CCP PETMR)
(http://www.ccpsynerbi.ac.uk/).

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*/

/*!
\file
\ingroup PET

\author Evgueni Ovtchinnikov
\author Richard Brown
\author SyneRBI
*/
#include <iostream>
#include <cstdlib>

#include "stir/common.h"
#include "stir/Verbosity.h"

#include "sirf/STIR/stir_x.h"
#include "sirf/common/getenv.h"
#include "sirf/common/iequals.h"

using namespace stir;
using namespace sirf;

int test_SPECT()
{
  std::cout << "running test4.cpp...\n";
  try {

    std::string SIRF_path = sirf::getenv("SIRF_PATH");
    if (SIRF_path.length() < 1) {
      std::cout << "SIRF_PATH not defined, cannot find data" << std::endl;
      return 1;
    }

    float * ptr = 0;
    //crash
    *ptr = 1.F;

    std::string path = SIRF_path + "/data/examples/PET/mMR/";
    fix_path_separator(path);

    // Construct STIRImageData from VoxelsOnCartesianGrid
    Coord3DI image_size = {31, 111, 111};
    Coord3DF voxel_size = {3.375, 3, 3};
    IndexRange3D index_range(0, image_size.z() - 1,
                             -(image_size.y() / 2), -(image_size.y() / 2) + image_size.y() - 1,
                             -(image_size.x() / 2), -(image_size.x() / 2) + image_size.x() - 1);
    Coord3DF offset = {0.f, 0.f, 0.f};

    shared_ptr<Voxels3DF> im_sptr(new Voxels3DF(
                                                index_range,
                                                offset,
                                                voxel_size));
    im_sptr->fill(0.0);
    STIRImageData stir_im(im_sptr);

    // Load mMR sinogram
    const std::string f_mMR_template = path + "mMR_template_span11.hs";
    PETAcquisitionDataInFile mMR_template(f_mMR_template.c_str());
    std::shared_ptr<PETAcquisitionDataInMemory> acq_data_mMR_sptr(
                                                                  new PETAcquisitionDataInMemory(mMR_template));
    acq_data_mMR_sptr->fill(1.f);

    // Create mMR image
    BasicCoordinate<3, int> min_image_indices(make_coordinate(0,  -160, -160));
    BasicCoordinate<3, int> max_image_indices(make_coordinate(126, 159,  159));
    IndexRange<3> range = IndexRange<3>(min_image_indices,max_image_indices);
    std::shared_ptr<STIRImageData> im_mMR_sptr(
                                               new STIRImageData(Voxels3DF(
                                                                           acq_data_mMR_sptr->get_exam_info_sptr()->create_shared_clone(),
                                                                           range,
                                                                           CartesianCoordinate3D<float>(0.f,0.f,0.f),
                                                                           CartesianCoordinate3D<float>(2.03125f, 2.08626f, 2.08626f))));

    im_mMR_sptr->fill(1.f);
    PETAcquisitionModelUsingRayTracingMatrix acq_model;
    stir::Verbosity::set(0);
    std::cout << "\nSetting up NiftyPET acquisition model...\n";
    acq_model.set_up(acq_data_mMR_sptr, im_mMR_sptr);
    std::cout << "\nForward projecting with NiftyPET acquisition model...\n";

    std::shared_ptr<PETAcquisitionData> prj_sptr = acq_model.forward(*im_mMR_sptr);
    std::cout << "\nBack projecting with NiftyPET acquisition model...\n";
    im_mMR_sptr = acq_model.backward(*prj_sptr);
    std::cout << "\nNiftyPET test succeeded.\n";


    std::cout << "done with test4.cpp...\n";
    return 0;
  }
  catch (...)
    {
      return 1;
    }
}

int main()
{
	const int failed = test_SPECT();
	std::cout << failed << " tests failed\n";
        
	return failed==0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
