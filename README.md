# M51 Merger Simulation: Velocity Analysis 
Adaptable code to run a basic simulation, written originally by [José Flores Velázquez](https://jaf12.github.io/joseflores.github.io/) and adapted to fit my science question: 
how does initial velocity of the incoming galaxy NGC 5195 impact the structure and evolution of M51?

## Getting Started 
These are just suggestions based on what I've seen work on my machine. There is probably a better way to do things--I'm just 
sharing what has worked for me! 

### About the code files
`Velocity_Merger.py` should run without any installations and give you several folders with potential merger movies in .png files. 
The same goes for `Isolated_M51.py`, which is a pared down version of the former file. It produces .png files illustrating an 
isolated evolution of M51, without the merger interaction. This is for comparison to the merger simulations and for creating 
a perturbation metric. This metric is visualized through the file `analysis.py`. This is to be run *_after_* the merger and 
isolation files have already run -- it requires positional data calculated within these files. `analysis.py` has two output 
options: a singular plot that illustrates the changing perturbation ratio over time, or a specified number of files based upon
the length of the simulation to create a .mp4 file that can be stitched with a merger video. 


### Using FFmpeg and creating videos from the code output
`Velocity_Merger.py` should run without any installations and give you several folders with potential merger movies. 
In order to turn the output .pdf files into the merger .mp4 file, you need to download [FFmpeg](https://www.ffmpeg.org/). 
Once FFmpeg is all downloaded, navigate into one of the directories with the .pdf files, and type the following: <br />

`ffmpeg -i frame%d.png output.mp4`

This will create a movie of the frames (in numerical order) with the name `output.mp4`. Now, if you find that you're having trouble with
the quality of the .mp4 output file, try running this to make the file instead: <br /> 

`ffmpeg -r 25 -i frame%d.png -vb 20M output.mpg`

If you're interested in stitching two simulation .mp4 files together, you first need to convert the files to the same size. 
You can do this through QuickTime (if you have MacOS) File > Export as > 1080p... and save the desired files to be 
stitched this way. You can also do this using the following command with [`avconv`](https://libav.org/avconv.html) (`brew install libav` first): 

`avconv -i input.mp4 -s hd1080 -c:v libx264 output.mp4`

Either of these should ensure both .mp4 files to be stitched together are the same size. 

## Authors

- Linnea Dahmen - adaptive and analysis work - Pomona College 
- José Flores Velázquez - initial merger code - UC Irvine 

## Acknowledgements 

- Professor Jorge Moreno - Simulation expertise and advising - Pomona College 
- Professor Gordon Stecklein - Project advising and research expertise - Pomona College 
- Professor Richard Mawhorter - Research expertise and support - Pomona College 
