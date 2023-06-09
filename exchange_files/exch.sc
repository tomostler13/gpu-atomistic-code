exchange:
{
    //This is the exchange interaction with the largest number of shells
    MaxShells=1;
};
exchange_0_0:
{
    Num_Shells = 1;
    //Units of energy
    units = "joules";
    Shell1Vec = [ 0.5 , 0.0 , 0.0 ];
    NumInt1 = 6;
    //Exchange Tensor First number is shell and second is row of tensor
    J1_1 = [ 5e-21 , 0.000 , 0.000 ];
    J1_2 = [ 0.000 , 5e-21 , 0.000 ];
    J1_3 = [ 0.000 , 0.000 , 5e-21 ];
};
exchange_0_1:
{
    Num_Shells = 1;
    //Units of energy
    units = "joules";
    Shell1Vec = [ 0.0 , 0.5 , 0.0 ];
    NumInt1 = 6;
    //Exchange Tensor First number is shell and second is row of tensor
    J1_1 = [ -1.09e-21 , 0.000 , 0.000 ];
    J1_2 = [ 0.000 , -1.09e-21 , 0.000 ];
    J1_3 = [ 0.000 , 0.000 , -1.09e-21 ];
};
exchange_1_1:
{
    Num_Shells = 1;
    //Units of energy
    units = "joules";
    Shell1Vec = [ 0.0 , 0.5 , 0.0 ];
    NumInt1 = 6;
    //Exchange Tensor First number is shell and second is row of tensor
    J1_1 = [ 1.26e-21 , 0.000 , 0.000 ];
    J1_2 = [ 0.000 , 1.26e-21 , 0.000 ];
    J1_3 = [ 0.000 , 0.000 , 1.26e-21 ];
};
exchange_1_0:
{
    Num_Shells = 1;
    //Units of energy
    units = "joules";
    Shell1Vec = [ 0.0 , 0.5 , 0.0 ];
    NumInt1 = 6;
    //Exchange Tensor First number is shell and second is row of tensor
    J1_1 = [ -1.09e-21 , 0.000 , 0.000 ];
    J1_2 = [ 0.000 , -1.09e-21 , 0.000 ];
    J1_3 = [ 0.000 , 0.000 , -1.09e-21 ];
};
