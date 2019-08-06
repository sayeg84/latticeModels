cd "src" 
for filename in *.jl; do
     echo "testing $filename" 
     julia -q $filename > /dev/null
done
for filename in *.py; do
     echo "testing $filename" 
     python $filename
done
cd ".."