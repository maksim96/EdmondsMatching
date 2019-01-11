for file in graphs/*.dmx
do
   echo $file
  ./EdmondsMatching $file
done
