The basic structure of the input form for calculations

```mermaid
   graph TD;
       A[Experiment]==>B[Sample]==>C(((Hamiltonian / Input Structure)))==>D[Sites];
       C==>E[Couplings]-.->D;
       A-->F(Magnetic Field);
       A==>G[Instrument]-->H(Resolution);
       B-->Ba[Sample Type]-->Bb(Q Points);
       
```

If we want to perform a fit to data rather than just a calculation, we can plug this into the Experiment class:

```mermaid 
graph TD;
    A[Experiment]==>B{Fit};
    C[Data]==>B;
    C-->D(Experimental Q); 
```

where the fit is over chosen parameters in the sample.

The input structure in the circle is generated from a more detailed, dynamic, user-friendly structure.
The MutableStructure class is responsible for making the structure user friendly. It is used internally by MATLAB style interface and provides the model for the GUI. 
The solid lines show the structure, and the dotted lines show actions that need to be taken to keep the stucture synchronised and valid.

```mermaid
  graph TD;
      A1@{shape: curv-trap, label: "MATLAB Style Interaction"}==>B[MutableStructure];
      A2@{shape: curv-trap, label: "GUI"}==>B
      B==>C[Sites]==>Ca[Independent Sites];
      C==>Cb[Implied Sites];
      B==>D[Symmetry];
      D==>Da[UnitCell];
      D==>Db[Supercell];
      D==>Dc[Crystal Type]==>Dd[Bravias]==>De[Space Group]==>Df[Magnetic Group];
      B==>E[Coupling Groups]-->F(((Generate Input Structure)));
      C-->F;
      D-.->G{Symmetry Check};
      C-.->G;
      C -.-> E;
      D -.-> C;
      Da -.-> C;
      Db -.-> C;
      Dc -.-> C;
      Dd -.-> C;
      De -.-> C;
      Df -.-> C;
      
```
