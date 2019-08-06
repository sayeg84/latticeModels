if ! hash julia
then
    echo "ERROR: Julia executable not found in PATH"
    exit 1
fi
if ! hash python
then
    echo "ERROR: Python executable not found in PATH"
    exit 1
fi
echo ""
echo ""
echo "Installing Julia dependencies ... "
echo ""
echo ""
julia -e "using Pkg ; Pkg.instantiate() "
echo ""
echo ""
echo "Done "
echo ""
echo ""
echo ""
echo ""
echo "Installing Python dependencies ... "
echo ""
echo ""
pip install -r REQUIREpython
echo ""
echo ""
echo "Done "
echo ""
echo ""
