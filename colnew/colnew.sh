#!/bin/bash
#
mkdir temp
cd temp
~/binc/f77split ../colnew.f
#
for FILE in `ls -1 *.f`;
do
  gfortran -c -w $FILE
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
done
rm *.f
#
ar qc libcolnew.a *.o
rm *.o
#
mv libcolnew.a ~/libf77
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/libcolnew.a."
