# Relationship-of-EMG-Signal-and-Muscular-Force
This is my group’s submission for the project of the course BME3121 Biomedical Signals and Systems in City University of Hong Kong. 
    
## Analysis of the Relationship Between Parameters of the EMG Signal and Muscular Force
The main objective of this project is to analyze the relationship between parameters of the EMG signal and muscular force by developing a method for automatic identification of portions (segments) corresponding to each level of contraction within which the force remains close to the corresponding peak values. 
    
Therefore, we created an algorithm based on finding local maxima and local minima to segment the signals. This successfully segmented the force signals that contributed the most to the EMG signal. When implemented in force signals with abrupt slope change from contraction to relaxation, it has a good performance, corroborated by various parameters such as its mean frequency, median frequency, and Higuchi Fractal Dimension.
    
Description of Files: 
1. BME3121_Project.m: Main Code
2. EMGforce1.txt, EMGforce2.txt, EM_EMG_SQUEEZE1.txt, and EM_EMG_SQUEEZE2.txt: Example data (with 3 columns representing time, force, EMG signal respectively)
