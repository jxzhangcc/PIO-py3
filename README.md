# PIO-py3
**Principal Interacting Orbital (for python 3)**

## Requirement
- Python 3
- Numpy 1.18.5
- Gaussian 09 or 16
- NBO 6.0 or above (recommended but not necessary)

This combination has been well tested. Other versions are not guaranteed to work but are welcomed to test.

## Tutorial
The overall workflow of a PIO analysis consists of three parts:
- Run ordinary quantum chemistry calculation, i.e. computing the electronic structure of a molecule and optimizing its geometry if necessary
- Perform natural population analysis (NPA)
- Perform principal interacting orbital (PIO) analysis

For Gaussian users, there is a built-in NBO 3.0 program, thus an individual NBO program is not necessary for performing PIO analysis. However, NBO 6.0 or above is usually believed to have superior performance compared to NBO 3.0 in terms of NPA charge, especially for transition metal centers, thus is recommended.
Due to the lack of universal interface to NBO program and visualization tools, incorpration with other programs is not implemented yet.

1. Run Gaussian calculation with NBO analysis

    **Input**

    - *CH4.gjf*
    ```
    %chk=CH4.chk
    #p opt freq b3lyp/6-31G* pop=nboread
    
    CH4
    
    0 1
     C                 -0.00000000    0.00000000    0.00000000
     H                  0.00000000    0.00000000    1.07000000
     H                 -0.00000000   -1.00880567   -0.35666667
     H                 -0.87365134    0.50440284   -0.35666667
     H                  0.87365134    0.50440284   -0.35666667
    
     $nbo FILENAME=CH4 AONAO=W33 FNAO=W61 DMNAO=W71 $end
     ```
    Shown above is a typical Gaussian input file for PIO analysis. The most important keyword is `pop=nboread` along with a separated line in the end of the file `$nbo FILENAME=CH4 AONAO=W33 FNAO=W61 DMNAO=W71 $end` which specifies the keywords for NBO program to save necessary infomation in local storage.

    **Output**
    - *CH4.chk*

        Gaussian checkpoint file

    - *CH4.33*

        NBO temp file storing NAO coefficients

    - *CH4.61*

        NBO temp file storing NAO density matrix

    - *CH4.71*

        NBO temp file storing NAO Fock matrix

     **Notes**
     - `FILENAME=CH4` specifies the filename of 33/61/71 files, otherwise the default filenames are FILE.33, etc.`
     - Fock matrix is not always available from Gaussian depending on different calculation objectives, in which case the subsequent procedure should still work except that PIO energies are no longer available in the final output.
     - For users who own and prefer to use NBO 6.0 or newer versions, `pop=nbo6read` can be used to call NBO 6.0 program. Alternatively one might use the NBO keyword `archive` to generate a `FILE.47` file and run external NBO program separately. In the latter case, be reminded to add `AONAO=W33 FNAO=W61 DMNAO=W71` in the NBO keyword line of the 47 file.

2. Prepare input files and perform PIO analysis
    Run the following commands in a shell environment.
    ```
    formchk CH4.chk CH4.fchk
    cat CH4.33 CH4.61 CH4.71 > CH4.49
    python PIO.py CH4.fchk
    ```
    These commands do three things:
    - Convert Gaussian checkpoint file to formcheck file for ease of visualization
    - Concatenate the NBO output files into one single file. Note that this file MUST have the same filename as the formcheck file so that the PIO program can properly find it.
    - Do PIO analysis

    In the last step, the PIO program will then request for a fragmentation as input to specify two groups of atoms with following prompt:

    ```
    $ Please input the atom ID of two fragments: (e.g. 1-5,8,13 6-7,9-12)
    $
    ```
    Two groups of atom IDs should be input here separated by a space. Numbers in each group are separated by a comma. Hyphen is supported for sequential numbers. Atom numbering starts from 1. Complete fragmentation is always recommended (i.e. the specified two groups cover all the atoms present in the system). Incomplete fragmentation will lead to absence of mathematical elegance but is still meaningful if you really want to do it.

    **Output**
    - *CH4_pio.txt*

        PIO log file, containing basic information of the PIOs of the system subject to the input fragmentation

    - *CH4_pio.fchk*

        Gaussian FormCheck file containing PIOs labeled as in the txt file, could be visualized by GaussView and other compatable orbital visualization softwares

    - *CH4_pimo.fchk*

        A similar Gaussian FormCheck file containing PIMOs whose ordering is same as that of PIOs

## Related publication
Original method of PIO: doi.org/10.1002/chem.201801220

Extension to spin-polarized systems: doi.org/10.1039/D0CP00127A

A recent review on PIO: doi.org/10.1002/wcms.1469

## Possible issues
- Assertion error:
```
  File "/home/user/softwares/PIO-py3/PIO/parse49.py", line 152, in find_n_from_n_tri
    assert M == i * (i + 1) // 2
AssertionError
```
This error arises because PIO program has not read enough data from files. One should first check if NBO analysis has normally terminates.
There are however some rare cases, usually when the basis set is very large (e.g. for large systems or with very diffuse orbitals), in which case the matrix files (33/61/71) are not properly printed with desired format, such as the following case:
```
  -136.863913355   -5.050851215  -94.820013190  115.472313795 -269.393527639
   -26.244735693 -144.731687559  160.791510377   45.250862633-1373.922308535
```
In NBO program, the matrix output has a fixed column width which is not enough for numbers beyond 1000. In this case, there will be no space between adjacent matrix elements and consequently PIO program cannot properly identify them as two numbers. To solve this issue, just make up the missing space. The PIO program simply uses space as delimiter so column width is not a concern.

- Key error:
```
  File "/home/user/softwares/PIO-py3/PIO/PIO.py", line 48, in PIO
    naoao = raw['NAOAO']
KeyError: 'NAOAO'
```
This error arises because PIO program didn't find the AO to NAO transformation matrix in the .49 file. Again, one should first check if NBO analysis has normally terminates.
However, incompatible Gaussian and NBO versions can also lead to this situation. To check if this is the case, open the .33 file and see if the file header is missing:
```
# NBO4
# NAOs in the AO basis:
# -------------------------------------------------------------------------------
    -0.032042813    1.024459509    0.022315269    0.001689888   -0.001077316
    ...
```
The lines starting with # should be present in most cases. But based on our tests, NBO3.0 built in G09 will NOT output this header. As NBO program uses this header to recognize which matrix is which, the PIO program cannot proceed.

In this case, you need to manually add back these three lines to the beginning of the .33 file, regenerate .49 file and rerun PIO program.
