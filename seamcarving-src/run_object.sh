g++ seam_carving_object.cpp -o seam_carving_object -Wall -O2 -std=c++11 -DSHOWSEAM

./seam_carving_object ../Images/obj1.jpg 1 --protect ../Images/obj1_green.jpg --remove ../Images/obj1_red.jpg 
mv output.ppm output1.ppm
mv output_showseam.ppm output1_showseam.ppm

./seam_carving_object ../Images/obj2.jpg 1 --protect ../Images/obj1_green.jpg --remove ../Images/obj1_red.jpg 
mv output.ppm output2.ppm
mv output_showseam.ppm output2_showseam.ppm
