ReadMe for NeEstimator Version 2.01 Minor Update May 2014

__________________________________

About the Minor Update of 2.01 May 2014

There was a typo in the v2.01 code that caused the lower jacknife confidence interval to be printed as ‘infinite’ when the upper interval was infinite. This has been corrected.

Another minor issue has been fixed also. Previously the software did not correctly print out information about monomorphic loci for temporal samples. This problem occurred when one population consisted of several samples and there were monomorphic loci present. This was not an issue for the single-sample methods (LD, Hets Excess and Molecular Coancestry).

____________________________________

You have downloaded the following files as a .zip package.

1)  NeEstimator 2.01.jar (the graphical user interface or GUI)
2)  Ne2.exe (executable for Windows), Ne2M (executable for mac) and Ne2L (executable for Linux).
    Users need to implement the correct version for their operating system.
3)  8Ne50.dat and 8Ne50.gen (test input files in two formats)
4)  8Ne50Ne.txt (example of a basic output file)
5)  Help 2.01 folder containing a help file as a .pdf and .html. Place the help
    folder in the same directory (folder) as the .jar and executables and the
    help file will launch on command from within the GUI in user's web browser.
6)  Accessory files (Mathematical details of methods implemented by software
    (NeCalcul.pdf), and files for batch processing of input files
    (multi.txt, multiplus.txt, and common.txt).

____________________________________

To run the software :-

1. Unpack (expand) the zip package.
2. Place into a suitable folder (directory) on your hard drive.
3. Run the NeEstimator V2 software by starting the graphical user interface as follows:
   * Windows or Mac Users: Double click on the NeEstimator.jar program.
     (There may be a lag on the Mac before the GUI appears)
   * Linux Users: From the command line execute:  "java -jar ./NeEstimator.jar". 
4. Load an input file using buttons on the interface. 
5. Click "Run NeEstimator".
6. On completion, look for output in the same directory as the input file. 
7. Click the "?" and "info" buttons within the user interface for details
   regarding numerous other options.
8. Read the help file for further information.

_____________________________________

About Mac OS and X11If running the Mac OS 10.8 or 10.9, you will need to download the latest version of XQuartz. Start up the .jar file (see above) and follow the on-screen OSX prompts to install XQuartz._____________________________________
NeEstimator development team
May 2014
