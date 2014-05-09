#!/bin/sh
for i in *.c *.h
do
echo $i
diff $i $1/$i
done
for i in */*.c */*.h
do
echo $i
diff $i $1/$i
done
