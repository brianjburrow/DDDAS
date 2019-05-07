Author: Brian Burrows
Email:  Brianjburrows@gmail.com

- This folder contains a small portion of (functioning) code that I wrote as a graduate student at Texas A&M.
- Most of the stuff that I tried was either not generally useful or not useful at all, so I have left those items out.
- To get the codes running, navigate to the parent folder containing the DDDAS folder.  Right click DDDAS in the MATLAB directory tree, and 
- click "add to path -> all selected files and folders".  Then all codes in Applications can be run in the usual way, and some functions from 

- The folders are:
-- Applications:
-- Methods:
-- oddballScripts:
-- Development:
-- utilityFunctions:

-- Applications Folder
- This folder contains direct applications of the Methods from the "Methods" folder, including the code that I used to generate the figures in my journal papers.
- Many of these are confusing, but they generate the plots from my papers.  This also contains a set of functions that I used to generate animations and visualitation
- tools for presentations.
- All of these files are directly executable, and rely on techniques stored in the subsequent folders.

-- Methods Folder
- This folder contains working packages for a few techniques from the literature, including 
-- 1. Constructing Transport Maps from a set of samples
-- 2. Bootstrap Particle filters and L2O Filters
-- 3. High Dimensional Model Representation codes for cut-HDMR and rs-HDMR (used for surrogate modeling aka emulation aka lots of things)
-- 4. Active Subspaces, which used for dimension reduction.
- Each method will contain test functions that show how to use the class object.
- Each method will contain pdf versions of the papers that I read when implementing those files, so that any weird variables can be tracked down.

-- oddballScripts Folder
- This is typically where I first write my visualization codes, and if I think they are useful, they get relocated.

-- Development:
- This is where I test new research ideas, new methods from journals that I read, etc.,
- Once the codes are functioning, they either get moved to the Applications folder (if I write a paper), the methods folder (if they are useful methods), or the utilityFunctions (if they are useful functions)
- Any codes in here should also have .pdf files of any necessary or useful papers.

-- utilityFunctions:
- these are typically weird numerical tricks (e.g., log-pdf functions), basis functions (i.e., multi-variate polynomial expansions), etc., that are required to get the above methods working.
- sometimes these are very weird plot functions that I wrote in order to make the above visualizations work.  I separate the necessary scripts from the executables so that it is easier to navigate the folders.
- Contains datasets, saved machine learning functions, etc., that are necessary to run some of the executables in "Applications".
- Contains toy problems for filtering, surrogate model testing, etc., with necessary .pdf files of the papers I pulled them from (or internet sources) when applicable