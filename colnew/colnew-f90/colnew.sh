#!/bin/bash
#
mkdir temp
cd temp
~/binc/f90split ../colnew.f90

#
for FILE in `ls -1 *.f90`;
do
  gfortran -c -ffree-form  $FILE
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
done
rm *.f90
#
ar qc libcolnew.a *.o
rm *.o
#
mv libcolnew.a ~/libf90
cd ..
rmdir temp
#
echo "Library installed as ~/libf90/libcolnew.a."
