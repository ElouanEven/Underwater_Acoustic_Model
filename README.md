# Modeling BELLHOP equation

Underwater acoustic using BELLHOP equation

Model include :
- Sound Speed Profile
- Reflexion for flat surface and bottom
- Varying sea bottom
- Absorption depending on frequency
- Intensity and transfer loss
- Receiver
    * detecting rays crossing
    * recording time and intensity delay
    * ploting only rays detect by receiver
- Hydrophone Array
- Segmentation in function
- Creating time signal at the source
- Add multiple receivers
- Recompose time signal on receivers and add noise
- Reduction of the range between source and receivers for time computation
- DaS Beamforming :
     * For direct rays
     * For reflected rays
- Optimization of the computational time by decorelating dt and Fs
- Return of the rays to search the position of the source
 
TO DO :
- Find the best spot of the source with time correlation
- Implement a threshold for max transfer loss
- Surface and bottom loss (tested but not clear)
- Use different input signal and their properties (sweep)
- Beamforming
    * MVDR
    * LS
    * MUSIC
- Optimize code
- Discover 2D*N and 3D reprensentation
- Try different form of arrays
- Applying spherical harmonic on hydrophone

Corrections :
- Pic curvature not C1 continuous 
- Scale everything down for time computation + Modification of variables physical nonsense


