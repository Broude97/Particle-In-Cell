// 1 VARIABLES
// 1.1 DEFAULT VALUES
/* Default values can be changed according to experimental
 conditions and
 * user's preferences */
  xy_scale = 200;		// Pixel size (nm/px)
  thd_lower = 650;		// Lower threshold to segment the particles
  obj_min_area = 10;		// Minimum size of particles (voxels)
  obj_max_area = 0;		/* Maximum size of particles (voxels)
  						set it to 0 (zero) for unlimited size:
  						0 => obj_max_vol = x * y * nSlices */
  mean_id = 1;			/* Mean intensity of a single particle.
					   	It corresponds to the mean integrated
				   		density (px values) and can be measured
				   		via Routine 4 - Calibration */
  cell_type = "Cell";	// Cell type
  nps_type = "NP";		// Particle type
  w = 6.0;				// Width of the cell membrane region (px)
  w_offset = 0;			/* Off-set for the cell membrane region.
				   		Negative = inwards, Positive = Outwards
				   		Only availabe for membrane region-3D */

// 1.2 SOME OTHER IMPORTANT VARIABLES
  var cell_vol_memb;		//Missing Variable
  var area_u;			//Missing Variable
  var area;			// Surface area of a slice
  var area_basal;	// Surface area - basal membrane
  var area_diff;	// Area difference between two adjacent slices
  var area_ratio;	// Ratio between perimeter of a slice and area_diff
  var area_region;	// Area of target region - routines 4 and 5
  var area_apical;	// Surface area of the cell (except basal membrane)
  var area_top;		// XY Area of the top most slice of the cell
  var bit8;			// Boolean value to verify if the upload image is binary
  var bg;			// Pixel value for background subtraction
  var bg_mode;		/* Boolean value for using or not slice_mode_mean
  		  	  		 as a constant value for background subtraction */
  var bx;			// Array for objects' bounding box in X
  var by;			// Array for objects' bounding box in Y
  var bz;			// Array for objects' bounding box in Z
  var bxx;			// Array for objects' width in X
  var byy;			// Array for objects' height in Y
  var bzz;			// Array for objects' depth in Z
  var cell_basal;		// Slice - central position of the basal membrane
  var cell_basal_end;	// Slice - end of the basal membrane region
  var cell_basal_p;		// cell_basal provesional (before substack)
  var cell_mask;		// File name: image mask of the cell
  var cell_outlines;	// File name: image of the cell with outlines
  var cell_name;		// File name w/o extension: loaded image of the cell
  var cell_path;		// File path: loaded image of the cell
  var cell_sub;			// File name: substack - image of the cell
  var cell_top;			// Slice for the top of the cell
  var cell_top_p;		// cell_top provesional (before substack)
  var cell_area;			// Area of the cell (membrane region + inside)
  var cell_area_file;	// File name: image of the cellular volume
  var cell_area_filled;	// Used to calculate the volume inside the cell (due to offset)
  var cell_area_in;		// Area inside the cell
  var cell_area_in_u;
  var color_cell;			// Color of the cell - merged image
  var color_nps_in;			// Color of intracellular particles - merged image
  var color_nps_memb;		// Color of apical membrane particles - merged image
  var color_nps_memb_basal;	// Color of basal membrane particles - merged image
  var color_nps_out;		// Color of  extracellular particles
  var dir_res;		// Directory for results
  var imj_fiji;		// Boolean for using Fiji macro language
  var imj;			// Boolean for using ImageJ macro language
  var exp_name;		// Identification of the experiment
  var locator;		// = 1 if a particle is detected, = 0 if not
  var memb_2D;		// Boolean for using the membrane region - XY-plane
  var memb_3D;		/* Boolean for using the membrane region-3D -
XY, YZ- and ZX-planes */
  var nps_bottom;	// Last slice in which particles are detected
  var nps_bottom_p;	// nps_bottom provesional (before substack)
  var nps_in;		// Total number of particles inside the cell
  var nps_memb;		// Total number of particles in the apical mebrane region
  var nps_memb_basal;	// Total number of particles in the basal mebrane region
  var nps_name;		// File name w/o extension: loaded image of particles
  var nps_path;		// File path: loaded image of particles
  var nps_seg;		// File name: image mask of segmented particles
  var nps_smo;		// File name: image mask of smoothed particles
  var nps_sub;		// File name: substack - image of the particles
  var nps_top;		// First slice in which particles are detected
  var nps_top_p;	// nps_top provesional (before substack)
  var nps_total;	// Total number of particles
  var nps_uptake;	// File name: colored image of particles and cell
  var obj_c;		// Counter for objects with one or more particles
  var obj_c_in;		// Counter for intracellular objects
  var obj_c_memb;	// Counter for apical membrane associated objects
  var obj_c_memb_basal;	// Counter for basal membrane associated objects
  var obj_id;		// Array for objects' IntDens
  var obj_in;		// Array for the nr. of particles inside the cell
  var obj_memb;		// Array for the nr. of particles in the apical membrane region
  var obj_memb_basal;	// Array for the nr. of particles in the basal membrane region
  var obj_nr_nps;	// Array for the nr. of particles in each object
  var obj_nr;		// Number of objects in the 3DOC results table
  var obj_out;		// Array for the nr. of particles outside the cell
  var obj_vol;		// Array for objects' volume (in voxels)
  var obj_x;		// Array for objects' position X
  var obj_y;		// Array for objects' position Y
  var obj_z;		// Array for objects' position Z
  var rev;			// Boolean for reversing the stacks
  var seg;			// Segmentation strategy
  var slice_cal;	// Number of slices to be analyzed - Routine 4
  var slice_max;	// Slice of maximum IntDens
  var slice_mode;	// Array with the modal pixel values of every slices
  var slice_mode_mean;	// Mean pixel value of slice_mode
  var slice_nr;			// Number of slices of the loaded image stack
  var t_hour;			// Variable for time - report
  var t_minute;			// Variable for time - report
  var thd_c_higher;		// Final higher thresold for segmenting the cell
  var thd_c_lower;		// Final lower thresold for segmenting the cell
  var thd_higher		// Higher threshold for particles. Default = pow(2, bitDepth())
  var thd_m1_higher;	// Automatic higher thresold for segmenting the cell
  var thd_m1_lower;		// Automatic lower thresold for segmenting the cell
  var thd_sd; 			// Array with the stdev of background for every slices
  var thd_sd_mean;		// Mean pixel value of thd_sd
  var thd_set;			// Boolean value for using or not thr_lower
  var total_id;			// Total IntDens
  var total_id_in;		// Total IntDens - intracellular particles
  var total_id_memb;	// Total IntDens - apical membrane associated particles
  var total_id_memb_basal;	// Total IntDens - basal membrane associated particles
  var area_unit;			/* Area pixel unit (nm^2)
  			  				area_unit = xy_scale * xy_scale */
  var w_slices;		// Width of the membrane region in z direction
  var w1;			// Width of the frame

  var ImCounter;

  // Precise Quantification Parameters
  var FluorescenseRatio		//Ratio of fluorescence intensity of the fluorophore outside the cell over that inside the cell. 
  var CalibrationFactor		//Ratio of mean flurorescence between micro and nanoparticles measured by fluorimetry

 	/*For more information refer to 
 	Gottstein, C., Wu, G., Wong, B. J., & Zasadzinski, J. A. (2013). 
 	Precise quantification of nanoparticle internalization. 
 	ACS nano, 7(6), 4933-4945
 	*/
ImCounter = 0;
do {
	ImCounter = ImCounter + 1;
	// 2 INITIAL SETTINGS
	  requires("1.49p");
	  run("Close All");
	  if (isOpen("Log") == true) {
		selectWindow("Log");
		run("Close");
	  }
	  run("Clear Results");
	  run("Colors...", "foreground=white background=black selection=yellow");
	  run("Line Width...", "line=1");
	  setForegroundColor(255, 255, 255);
	  setBackgroundColor(0, 0, 0);
	  setOption("BlackBackground", true);
	  run("Options...", "iterations=1 count=1 black edm=Overwrite");
	  roiManager("Associate", "true");
	  selectWindow("ROI Manager");
	  run("Close");
	  ChangeParam = true;
	  if (ImCounter > 1) {
	  	Dialog.create("Parameters");
	  	Dialog.addCheckbox("Change Parameters?", false);
	  	Dialog.show();
	  	ChangeParam = Dialog.getCheckbox();
	  }
		
	  if (ImCounter == 1 || ChangeParam) {
		// 3 ROUTINE SELECTION
		  do {
			Dialog.create("Particle_in_Cell-3D");
			Dialog.setInsets(10,10,0);
			Dialog.addMessage("ROUTINE SELECTION");
			Dialog.setInsets(5,10,0);
			Dialog.addMessage("Particle Uptake by Single Cells");
			Dialog.setInsets(0,20,0);
			Dialog.addCheckbox("1. Qualitative - Visualization of the intra"
			+"cellular distribution of particles", false);
			Dialog.addCheckbox("2. Semi-Quantitative - Visualization and "
			+"intensity of particles", false);
			Dialog.addCheckbox("3. Quantitative - Visualization and "
			+"absolute number of particles", false);
			Dialog.setInsets(5,10,0);
			Dialog.addMessage("Particle Characterization");
			Dialog.addCheckbox("4. Calibration - Intensity of single particles", false);
			Dialog.addCheckbox("5. Only Particles - Intensity and absolute number of "
			+"particles", false);
			Dialog.addMessage("Additional Options");
			Dialog.addCheckbox("Multiple Cells", false);
			Dialog.addCheckbox("Precise Quantification", false);
			Dialog.addCheckbox("Multiple Images", false);
			Dialog.addCheckbox("Control Cells", false);
			Dialog.addCheckbox("Background Substraction", false);
			Dialog.show();
			R1 = Dialog.getCheckbox();
			R2 = Dialog.getCheckbox();
			R3 = Dialog.getCheckbox();
			R4 = Dialog.getCheckbox();
			R5 = Dialog.getCheckbox();
			MultCells = Dialog.getCheckbox();
			PreciseQuant = Dialog.getCheckbox();
			MultIm = Dialog.getCheckbox();
			ControlCells = Dialog.getCheckbox();
			BackSub = Dialog.getCheckbox();
			imj_fiji = true;
			choice = R1 + R2 + R3 + R4 + R5;
			title = "Particle_in_Cell-3D";
			if (choice >= 2) {
				msg = "Select only one routine!";
				waitForUser(title, msg);
			}
			if (choice == 0) {
				msg = "Select at least one routine!";
				waitForUser(title, msg);
			}
		  } while (choice >= 2 || choice == 0);
	  }
	// 4 SELECTION OF FILES
	  run("Close All");
	  nps_path = File.openDialog("Select the image of the PARTICLES");
	  nps_name = File.nameWithoutExtension;
	  if (R1 || R2 || R3 == true) {
		cell_path = File.openDialog("Select the image of the CELL");
		cell_name = File.nameWithoutExtension;
	  }

	  if (ControlCells) {
	  	ControlCellsPath = File.openDialog("Select Control Cells Image");
	  	ControlCellsName = File.nameWithoutExtension;
	  }
	  
	  dir_res = getDirectory("Select a directory for RESULTS");
	
	// 5 REGION SELECTION (XY-PLANE)
	  if (R1 || R2 || R3 == true) {
		open(cell_path);
		run("Select None");
		run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
		getImageInfo();
		bit8 = getInfo("BitsPerPixel");
		if (bit8 == 8); {
			run("16-bit");}
		if (isOpen("Log") == true) {
			selectWindow("Log");
			run("Close");
		}
	
		// 5.1 SUGGESTION FOR THE Z-POSITION OF THE BASAL MEMBRANE
		slice_nr = nSlices;
		run("Set Measurements...", "  integrated stack redirect=None "
		+"decimal=0");
		run("Measure");
		slice_id_max = getResult("RawIntDen", 0);
		selectWindow("Results");
		run("Close");
		run("Enhance Contrast...", "saturated=0.8 use");
	  }
	  if (R4 || R5 == true) {
	
		// 5.2 SUGGESTION FOR THE Z-POSITION OF DEPOSITED PARTICLES
		open(nps_path);
		run("Select None");
		run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
		getImageInfo();
		bit8 = getInfo("BitsPerPixel");
		if (bit8 == 8); {
			run("16-bit");}
		if (isOpen("Log") == true) {
			selectWindow("Log");
			run("Close");
		}
		slice_nr = nSlices;
		run("Set Measurements...", "  integrated stack redirect=None "
		+"decimal=0");
		run("Measure");
		slice_id_max = getResult("RawIntDen", 0);
		selectWindow("Results");
		run("Close");
		run("Enhance Contrast...", "saturated=0.8 use");
	  }
	  if (R1 || R2 || R3 == true) {
	
		// 5.3 TARGET CELL
		setTool("freehand");
		title = "Particle_in_Cell-3D";
		msg = "TARGET CELL"
		+"\n \nIf the target cell is in contact with other cell(s) "
		+"\nselect a region so as to segregate her."
		+"\nDo nothing if it is already alone."
		+"\n \nClick 'OK' to go on.";
		waitForUser(title, msg);
		sel_type = selectionType();
		if (sel_type == -1) {
			run("Select All"); }
		roiManager("Add");
		roiManager("Save", dir_res+ "TARGET_REGION-" +nps_name+ ".zip");
		run("Clear Outside", "stack");
		run("Enhance Contrast...", "saturated=0.8 use");
		run("Select None");
		run("Set Measurements...", "  integrated stack redirect=None "
		+"decimal=0");
		run("Measure");
		slice_id_max = getResult("RawIntDen", 0);
		selectWindow("Results");
		run("Close");
		saveAs("tiff", dir_res + "CELL-" + cell_name);
		close();
		open(""+ dir_res + "CELL-" + cell_name +".tif");
		cell_sub = File.name;
		close();
	  }
	  if (R4 || R5 == true) {
	
		// 5.4 TARGET PARTICLES
		setTool("rectangle");
		title = "Particle_in_Cell-3D";
		msg = "TARGET REGION (XY-PLANE)"
		+"\n \nDraw a region of interest. "
		+"\nDo nothing if it is the whole area."
		+"\n \nClick 'OK' to go on.";
		waitForUser(title, msg);
		sel_type = selectionType();
		if (sel_type == -1) {
			run("Select All"); }
		roiManager("Add");
		roiManager("Save", dir_res+ "TARGET_REGION-" +nps_name+ ".zip");
		run("Clear Outside", "stack");
		run("Enhance Contrast...", "saturated=0.8 use");
		run("Set Measurements...", "area redirect=None decimal=0");
		run("Measure");
		area_region = getResult("Area", 0);
		selectWindow("Results");
		run("Close");
		run("Select None");
		run("Set Measurements...", "  integrated stack redirect=None "
		+"decimal=0");
		run("Measure");
		slice_id_max = getResult("RawIntDen", 0);
		selectWindow("Results");
		run("Close");
	  }
	
	 // 6 BACKGROUND SUBTRACTION
	  if (R1 || R2 || R3 == true) {
		open(nps_path);
		run("Select None");
		run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
		getImageInfo();
		bit8 = getInfo("BitsPerPixel");
		if (bit8 == 8); {
			run("16-bit");}
		if (isOpen("Log") == true) {
			selectWindow("Log");
			run("Close");
		}
		roiManager("Select", 0);
		run("Clear Outside", "stack");
		run("Enhance Contrast...", "saturated=0.8 use");
	  }
	  run("Select None");
	  saveAs("tiff", dir_res + "PARTICLES-" + nps_name);
	  close();
	  open(""+ dir_res + "PARTICLES-" + nps_name+".tif");
	  nps_sub = File.name;
	  if (obj_max_area == 0) {
		x = getWidth();
		y = getHeight();
		obj_max_area = x * y;
	  }
	
	  // 6.1 MODAL PIXEL VALUE SUGGESTED FOR BACKGROUND SUBTRACTION
	  slice_mode = newArray(nSlices);
	  run("Set Measurements...", "  modal redirect=None decimal=0");
	  roiManager("Select", 0);
	  run("Measure");
	  slice_mode[0] = getResult("Mode", 0);
	  Array.getStatistics(slice_mode, min_md, max_md, mean_md, sd_md);
	  slice_mode_mean = round(mean_md);
	  selectWindow("Results");
	  run("Close");

	  if (ControlCells) {
	  	open(ControlCellsPath);
		run("Threshold...");
		setAutoThreshold("IJ_IsoData dark stack");
		title = "Particle_in_Cell-3D";
		msg = "IMAGE MASK OF THE Control Cells";
		waitForUser(title, msg);
		getThreshold(thd_c_lower, thd_c_higher);
		run("Convert to Mask", " background=Dark black");
		run("Select None");
		selectWindow("Threshold");
		run("Close");
		run("Fill Holes", "slice");
		i = 0;
		do {
			setTool(8);
			title = "Particle_in_Cell-3D";
			msg = "JUST A FEW CLICKS!"
			+"\n \nClick on the white mask of the cell to use as control";
			waitForUser(title, msg);
			roiManager("Add");
			roiManager("Select", i);
			run("Interpolate", "interval=3 smooth");
			roiManager("Update");
			Dialog.create("Particles in Cell");
			Dialog.addCheckbox("More Cells?", false);
			Dialog.show();
			MoreCellsCont = Dialog.getCheckbox();
			i = i + 1;
		} while (MoreCellsCont);
	  	run("Set Measurements...", "area mean redirect=None decimal=0");
		roiManager("deselect");
		roiManager("combine");
		open(ControlCellsPath);
		roiManager("Select",0);
		run("Measure");
		bg = getResult("Mean", 0);
		selectWindow("Results");
		run("Close");
		bg_mode = 1;
		roiManager("reset");
		run("Close All");
	  }
	  else {
	  bg_mode = getBoolean("BACKGROUND SUBTRACTION"
	  +"\n \nThe calculated image background is "+slice_mode_mean+"."
	  +"\n \nDo you want to subtract this value from the image?"
	  +"\n \n  'YES' to go on"
	  +"\n  'NO' to set another value"
	  +"\n  'CANCEL' to exit");
	  }
	  // 6.2 SUBTRACTION
	  if (bg_mode == 0) {
		Dialog.create("Particle_in_Cell-3D");
		Dialog.addNumber("Background to be subtracted:"
		+"", bg, 0, 6,"px value");
		Dialog.show();
		bg =  Dialog.getNumber();
	  }
	  if (bg_mode == 1) {
		bg = slice_mode_mean;}
	  open(""+ dir_res + "PARTICLES-" + nps_name+".tif");
	  selectWindow(nps_sub);
	  run("Select None");
	  run("Subtract...", "value=&bg stack");
	  run("Enhance Contrast...", "saturated=0.8 use");
	  run("Select None");
	  saveAs("tiff", dir_res + "PARTICLES-" + nps_name);
	  run("Close All");
	  roiManager("Reset");
	
	// 7 ANALYSIS PARAMETERS
	if (ImCounter == 1) {
		//7.1 Precise Quantification Parameters
	
			if (PreciseQuant) {
				Dialog.create("Precise Quantification Parameters");
				Dialog.addNumber("Ratio of fluorescence Intensity", 1);
				Dialog.addMessage("Ratio of fluorescence intensity of the fluorophore");
				Dialog.addMessage("\noutside the cell over that inside the cell.");
				Dialog.addNumber("Calibration factor", 1);
				Dialog.addMessage("Ratio of mean flurorescence between");
				Dialog.addMessage("\nmicro and nanoparticles measured by fluorimetry.");
				Dialog.addMessage("\nFor more information refer to") 
	 			Dialog.addMessage("\nGottstein, C., Wu, G., Wong, B. J., & Zasadzinski, J. A. (2013).");  
	 			Dialog.addMessage("\nPrecise quantification of nanoparticle internalization."); 
	 			Dialog.addMessage("\nACS nano, 7(6), 4933-4945");
	 			Dialog.show();
	 			FluorescenseRatio = Dialog.getNumber();
	 			CalibrationFactor = Dialog.getNumber();
			}
		}
	  thd_set = 0;
	  do {
	  	if (ImCounter == 1) {
			run("Close All");
			roiManager("Reset");
			Dialog.create("Particle_in_Cell-3D");
			Dialog.addMessage("Please Enter Values for...",12);
			Dialog.addMessage("\nIDENTIFICATION",12);
			if (R1 || R2 || R3 == true) {
				Dialog.addString("Cell type:", cell_type, 18);}
			Dialog.addString("Particle type:", nps_type, 18);
			Dialog.addString("Experiment:", nps_name, 18);
			Dialog.addMessage("\nANALYSIS PARAMETERS",12);
			if (R1 == false) {
				Dialog.addNumber("XY-scale:", xy_scale, 1, 6,"nm/px");
			}
			if (R1 || R2 || R3 == true) {
				Dialog.addNumber("Width - membrane region:", w, 1, 6,"px, XY-scale");
				Dialog.addMessage("Off-set positive=> outwards, negative=> inwards",12);
				Dialog.addMessage("")
				}
			Dialog.addNumber("Lower threshold:", thd_lower, 0, 15, "px value");
			Dialog.addNumber("Minimum area for objects:"
			+"", obj_min_area, 0, 15, "voxels");
			Dialog.addNumber("Maximum area for objects:"
			+"", obj_max_area, 0, 15, "voxels");
			if (R3 || R5 == true) {
				Dialog.addNumber("Mean intensity of single particles:"
				+"", mean_id, 0, 15, "px value");
			}
			Dialog.addCheckbox("Exclude objects on edges", false);
			Dialog.addMessage("\nCOLOR CODING",12);
			if (R1 || R2 || R3 == true) {
				Dialog.addChoice("Cell:"
				+"", newArray("Green", "Red", "Cyan", "Magenta", "Grays"));
				Dialog.addChoice("Particles inside:"
				+"", newArray("Magenta", "Green", "Red"));
				Dialog.addChoice("Particles apical membrane:"
				+"", newArray("Yellow", "Cyan"));
				//Dialog.addChoice("Particles Outside:"
				//+"", newArray("Do Not Show", "Grays"));
			}
			if (R4 || R5 == true) {
				Dialog.addChoice("Particles:"
				+"", newArray("Magenta", "Green", "Red"));
			}
			Dialog.addMessage("\n \nHINTS:",12);
			Dialog.addMessage("You can try out and come back to this window",12);
			Dialog.show();
			if (R1 || R2 || R3 == true) {
				cell_type = Dialog.getString();}
			nps_type = Dialog.getString();
			exp_name = Dialog.getString();
			if (R1 == false) {
				xy_scale =  Dialog.getNumber();
				area_unit = xy_scale * xy_scale;
			}
			if (R1 || R2 || R3 == true) {
				memb_2D = true;
				w =  Dialog.getNumber();
			}
			thd_lower = Dialog.getNumber();
			obj_min_area = Dialog.getNumber();
			obj_max_area = Dialog.getNumber();
			if (R3 || R5 == true) {
				mean_id = Dialog.getNumber();}
			if (R1 || R2 || R4 == true) {
				mean_id = 1;}
			exclude = Dialog.getCheckbox();
			if (R1 || R2 || R3 == true) {
				color_cell = Dialog.getChoice();
				color_nps_in = Dialog.getChoice();
				color_nps_memb = Dialog.getChoice();
				color_nps_out = "Do Not Show";
			}
			if (R4 || R5 == true) {
				color_nps_in = Dialog.getChoice();
				MM = 1;
			}
	  	}
		// 7.1 THRESHOLD FOR THE PARTICLES
		// 7.1.1 SEGMENTATION
		open(nps_sub);
		roiManager("Open", dir_res + "TARGET_REGION-" + nps_name + ".zip");
		roiManager("Select", 0);
		run("Gaussian Blur...", "sigma=1.0 stack");
		run("Select None");
		saveAs("tiff", dir_res + "PARTICLES_SEGMENTED-" + nps_name);
		close();
		open(""+ dir_res + "PARTICLES_SEGMENTED-" + nps_name +".tif");
		nps_seg = File.name;
		roiManager("Reset");
		thd_higher = pow(2, bitDepth());
		setThreshold(thd_lower, thd_higher);
		run("Create Selection");
		sel_type = selectionType();
		if (sel_type == -1) {
			run("Select All");}
		run("Enlarge...", "enlarge=3");
		roiManager("Add");
		roiManager("Save", dir_res + "OBJECTS_REGION-" + nps_name + ".zip");
		roiManager("Reset");
		run("Close All");
	
		// 7.1.2 SUGGESTION
		open(nps_sub);
		roiManager("Open", dir_res + "TARGET_REGION-" + nps_name + ".zip");
		roiManager("Select", 0);
		thd_sd = newArray(nSlices);
		run("Set Measurements...", "  standard redirect=None decimal=0");
		for (i = 0; i < nSlices; i++) {
			setSlice(i + 1);
			run("Measure");
		}
		for (i = 0; i < nSlices; i++) {
			thd_sd[i] = getResult("StdDev", i);}
		Array.getStatistics(thd_sd, min_sd, max_sd, mean_sd, sd_sd);
		thd_sd_mean = round(mean_sd);
		selectWindow("Results");
		run("Close");
		close();
	
		// 7.1.3 SETTING THE VALUE
		open(nps_seg);
		roiManager("Select", 0);
		run("Threshold...");
		setAutoThreshold("Default dark stack");
		resetThreshold;
		setThreshold(thd_lower, thd_higher);
		wait(3000);
		title = "Particle_in_Cell-3D";
		msg = "THRESHOLD"
		+"\n \nVerify if the objects of interest are below the threshold."
		+"\n \nHINTS: \n-The threshold should be set as low as possible, "
		+"\nbut high enough to allow object segmentation."
		+"\n-Use the same threshold when comparing results from the same "
		+"\nset of experiments.";
		waitForUser(title, msg);
		getThreshold(thd_lower, thd_higher);
		thd_higher = pow(2, bitDepth());
		thd_set = getBoolean("THRESHOLD = "+thd_lower+" ?"
		+"\n 'YES' to continue"
		+"\n 'NO' to get a suggestion and then set a new value"
		+"\n 'CANCEL' to exit");
		if (thd_set == 0) {
			thd_lower = 4 * thd_sd_mean;
			msg = "SUGGESTION FOR THE THRESHOLD = "+4 * thd_sd_mean+""
			+"\n (4 * StdDev of image background)";
			waitForUser(title, msg);
		}
	  } while (thd_set == 0);
	  roiManager("Select", 0);
	  run("Clear Outside", "stack");
	  run("Select None");
	  saveAs("tiff", dir_res + "PARTICLES_SEGMENTED-" + nps_name);
	  roiManager("Reset");
	
	  if (R4 || R5 == true) {
	
		// 8 SUBSTACK SELECTION (Z-AXIS) - CALIBRATION & ONLY PARTICLES
		title = "Particle_in_Cell-3D";
		msg = "SUBSTACK SELECTION"
		+"\n \nPlease check the substack to be analyzed. "
		+"You will be soon asked to enter:"
		+"\n \n- The position of the cover glass - Suggestion: slice no. "+slice_max+""
		+"\n- The number of slices to be analyzed "
		+"\nOR"
		+"\n- The first slice"
		+"\n- The last"
		+"\n \nClick 'OK' to continue.";
		waitForUser(title, msg);
		Dialog.create("Particle_in_Cell-3D");
		Dialog.addMessage("SUBSTACK SELECTION");
		Dialog.addMessage("Please Enter the Values for...");
		Dialog.addNumber("Position of the cover glass - Slice:", slice_max);
		Dialog.addNumber("Number of slices:", 0);
		Dialog.addMessage("OR");
		Dialog.addNumber("First slice:", 1);
		Dialog.addNumber("Last slice:", slice_nr);
		Dialog.show();
		slice_max =  Dialog.getNumber();
		slice_cal =  Dialog.getNumber();
		nps_top_p =  Dialog.getNumber();
		nps_bottom_p =  Dialog.getNumber();
	
		// 8.1 SUBSTACK - PARTICLES
		if (slice_cal != 0) {
			nps_top_p = round(slice_max - (slice_cal/2));
			nps_bottom_p = round(slice_max + (slice_cal/2)) - 1;
		}
		selectWindow(nps_seg);
		run("Select None");
		run("Make Substack...", "slices="+nps_top_p+"-"+nps_bottom_p+"");
		saveAs("tiff", dir_res + "PARTICLES_SEGMENTED-" + nps_name);
		close();
		open(nps_sub);
		run("Select None");
		run("Make Substack...", "slices="+nps_top_p+"-"+nps_bottom_p+"");
		saveAs("tiff", dir_res + "PARTICLES-" + nps_name);
		cell_vol = nSlices * area_region;
		cell_vol_u = (area_unit * cell_vol)/1E9;
	
		// 8.2 UPDATING SLICE NUMBER AFTER SUBSTACK
		nps_bottom = nps_bottom_p - nps_top_p + 1;
		nps_top = 1;
		slice_max = slice_max - nps_top_p + 1;
		if (slice_cal == 0) {
			run("Select None");
			run("Set Scale...", "distance=0 known=0 pixel=1 unit=pixel");
			slice_nr = nSlices;
			run("Set Measurements...", "  integrated stack redirect=None "
			+"decimal=0");
			slice_id_max = 0;
			for(i = 0; i < slice_nr - 1; i++) {
				setSlice(i + 1);
				run("Measure");
				slice_id = getResult("RawIntDen", i);
				if (slice_id > slice_id_max) {
					slice_id_max = slice_id;
					slice_max = getResult("Slice", i);
				}
			}
			selectWindow("Results");
			run("Close");
		}
		run("Close All");
	  }
	  
	  if (R1 || R2 || R3 == true) {
	
		// 9 THREE-DIMENSIONAL RECONSTRUCTION OF THE CELL - ROIs
		// 9.1  IMAGE MASK OF THE CELL
		open(cell_sub);
		roiManager("Open", dir_res + "TARGET_REGION-" + nps_name + ".zip");
		roiManager("Select", 0);
		run("Gaussian Blur...", "sigma=1.0 stack");
		run("Select None");
		run("Duplicate...", "title=Mask-1 duplicate range=[0]");
		selectWindow(cell_sub);
		run("Duplicate...", "title=Mask-2 duplicate range=[0]");
		run("Tile");
	
		// 9.2 AUTOMATIC THRESHOLD FOR THE CELL
		selectWindow("Mask-1");
		roiManager("Select", 0);
		run("Threshold...");
		setAutoThreshold("IJ_IsoData dark stack");
		getThreshold(thd_m1_lower, thd_m1_higher);
		run("Convert to Mask", "method=IJ_IsoData background=Dark black");
		run("Select None");
	
		// 9.3 USER-SET THRESHOLD FOR THE CELL
		selectWindow("Mask-2");
		roiManager("Select", 0);
		setAutoThreshold("IJ_IsoData dark");
		title = "Particle_in_Cell-3D";
		msg = "IMAGE MASK OF THE CELL"
		+"\n \nIf 'Mask-1' is not satisfactory, set the threshold in 'Mask-2' "
		+"\nthat best defines the outlines of the cell.";
		waitForUser(title, msg);
		selectWindow("Mask-2");
	
		getThreshold(thd_c_lower, thd_c_higher);
		run("Convert to Mask", " background=Dark black");
		run("Select None");
		roiManager("Reset");
		selectWindow("Threshold");
		run("Close");
	
		// 9.4 IMAGE MASK OF THE CELL
		run("Tile");
		title = "Particle_in_Cell-3D";
		msg = "THE BEST IMAGE MASK"
		+"\n \nVerify which image mask can better represent the cell."
		+"\n \nClick 'OK' to continue.";
		waitForUser(title, msg);
		do {
			Dialog.create("Particle_in_Cell-3D");
			Dialog.addMessage("IMAGE MASK OF THE CELL");
			Dialog.addMessage("Select one option:");
			Dialog.addCheckbox("Mask-1 (Automatic)", false);
			Dialog.addCheckbox("Mask-2 ", false);
			Dialog.show();
			M1 = Dialog.getCheckbox();
			M2 = Dialog.getCheckbox();
			M = M1 + M2;
		} while (M != 1);
		if (M1 == true) {
			thd_c_lower = thd_m1_lower;
			thd_c_higher = thd_m1_higher;
			selectWindow("Mask-2");
			close();
			selectWindow("Mask-1");
		}
		if (M2 == true) {
			selectWindow("Mask-1");
			close();
			selectWindow("Mask-2");
		}
		run("Select None");
		saveAs("tiff", dir_res + "CELL_MASK_ORIGINAL-" + cell_name);
		close();
		open(""+ dir_res + "CELL_MASK_ORIGINAL-" + cell_name +".tif");
		cell_mask = File.name;
		run("Tile");
		
		// 9.5 SUBSTACK SELECTION (Z-AXIS) - PARTICLE UPTAKE
		// 9.5.1 SELECTION OF IMPORTANT SLICES
		cell_top_p =  1;
		cell_basal_p = 1;
		nps_top_p =  1;
		nps_bottom_p =  1;
	
		// 9.5.3 UPDATING SLICE NUMBER AFTER SUBSTACK
		cell_top = 1;
		cell_basal = 1;
		nps_bottom = 1;
		nps_top = 1;
		slice_max = 1;
	
		if (MultCells == false) {
		
			// 9.6.1 SEED FOR TRACKING THE POSITION OF THE CELL
			selectWindow(cell_mask);
			setSlice(cell_top);
			setTool("point");
			setSlice(cell_top);
			title = "Particle_in_Cell-3D";
			msg = "JUST A CLICK!"
			+"\n \nClick on the white mask of slice no. "+cell_top+""
			+"\nIt will create a point selection. Then click 'OK'.";
			waitForUser(title, msg);
			run("Draw", "stack");
			getSelectionCoordinates(xSeed, ySeed);
			run("Select None");
		
			// 9.6.2 SEGMENTATION STRATEGY - S1 OUTER-ROI
			selectWindow(cell_mask);
			rename("S1_Mask");
			run("Duplicate...", "title=S2_Mask duplicate range=[0]");
			setPasteMode("Add");
			selectWindow("S1_Mask");
			for(i = cell_top; i <= nps_bottom; i++) {
				if (i < cell_basal) {
					setSlice(i);
					run("Select None");
					run("Close-", "slice");
					run("Fill Holes", "slice");
					doWand(xSeed[0], ySeed[0]);
					run("Clear Outside", "slice");
					roiManager("Add");
					run("Select All");
					//Closing the mask in case it touches the frame
					run("Enlarge...", "enlarge=-4");
					run("Make Inverse");
					run("Copy");
					setSlice(i + 1);
					run("Paste");
				}
				if (i == cell_basal) {
					run("Select None");
					run("Close-", "slice");
					run("Fill Holes", "slice");
					doWand(xSeed[0], ySeed[0]);
					run("Clear Outside", "slice");
					roiManager("Add");
				}
				if (i > cell_basal) {
					j = cell_basal - cell_top;
					roiManager("Select", j);
					setSlice(i);
					roiManager("Add");
				}
			}
		}
	
		if (MultCells == true) {
			selectWindow(cell_mask);
			rename("S1_Mask");
			selectWindow("S1_Mask");
			run("Fill Holes", "slice");
			MoreCells = true;
			i = 0;
			run("Set Measurements...", "bounding redirect=None decimal=0");
			Roix = newArray;
			Roiy = newArray;
			RoiWidth = newArray;
			RoiHeight = newArray;
			setOption("ExpandableArrays", true);
			do {
				setTool(8);
				title = "Particle_in_Cell-3D";
				msg = "JUST A FEW CLICKS!"
				+"\n \nClick on the white mask of the cell of interest";
				waitForUser(title, msg);
				run("Measure");
				Roix[i] = getResult("BX", i);
				Roiy[i] = getResult("BY", i); 
				RoiWidth[i] = getResult("Width", i);
				RoiHeight[i] = getResult("Height", i);
				roiManager("Add");
				roiManager("Select", i);
				run("Interpolate", "interval=3 smooth");
				roiManager("Update");
				Dialog.create("Particles in Cell");
				Dialog.addCheckbox("More Cells?", false);
				Dialog.show();
				MoreCells = Dialog.getCheckbox();
				i = i + 1;
			} while (MoreCells);
			roiManager("deselect");
			roiManager("combine");
			run("Clear Outside", "slice");
			roiManager("add");
			roiManager("Save", dir_res + "S1-" + cell_name + ".zip");
			roiManager("Reset");
		}
	
		if (MultCells == false) {
			for(i = 0; i <= (nps_bottom - cell_top); i++) {
				roiManager("Select", i);
				run("Interpolate", "interval=3 smooth");
				roiManager("Update");
			}
	
			run("Select None");
			setSlice(slice_max);
			roiManager("Deselect");
			roiManager("Save", dir_res + "S1-" + cell_name + ".zip");
			roiManager("Reset");
		}
		
	
	
	
		// 9.6.3 SEGMENTATION STRATEGY - S2 OUTER-ROI
	/*	selectWindow("S2_Mask");
		for(i = cell_top; i <= nps_bottom; i++) {
			if (i <= cell_basal) {
				setSlice(i);
				run("Select None");
				run("Close-", "slice");
				run("Fill Holes", "slice");
				doWand(xSeed[0], ySeed[0]);
				run("Clear Outside", "slice");
				roiManager("Add");
				if (i != cell_basal) {
					run("Copy");
					setSlice(i + 1);
					run("Paste");
				}
			}
			if (i > cell_basal) {
				j = cell_basal - cell_top;
				roiManager("Select", j);
				setSlice(i);
				roiManager("Add");
			}
		}
		for(i = 0; i <= (nps_bottom - cell_top); i++) {
			roiManager("Select", i);
			run("Interpolate", "interval=3 smooth");
			roiManager("Update");
		}
		run("Select None");
		setSlice(slice_max);
		roiManager("Deselect");
		roiManager("Save", dir_res + "S2-" + cell_name + ".zip");
		roiManager("Reset");
	*/
		// 9.6.4 OUTLINES FOR S1 AND S2
		run("Close All");
		setForegroundColor(255, 255, 0);
		open(cell_sub);
		rename("S1-Outlines");
		roiManager("Open", dir_res + "TARGET_REGION-" + nps_name + ".zip");
		roiManager("Select", 0);
		run("Enhance Contrast...", "saturated=0.8 process_all use");
		run(color_cell);
		run("RGB Color");
		run("Select None");
		run("Duplicate...", "title=S2-Outlines duplicate range=[0]");
		roiManager("Reset");
		roiManager("Open", dir_res + "S1-" + cell_name + ".zip");
		selectWindow("S1-Outlines");
		roiManager("Select", roiManager("count") - 1);
		run("Draw", "slice");
		run("Select None");
		setSlice(slice_max);
		roiManager("Reset");
		/*
		roiManager("Open", dir_res + "S2-" + cell_name + ".zip");
		selectWindow("S2-Outlines");
		for(i = 0; i <= (nps_bottom - cell_top); i++) {
			roiManager("Select", i);
			run("Draw", "slice");
		}
		run("Select None");
		setSlice(slice_max);
		roiManager("Reset");
		run("Tile");
		*/
		// 9.6.5 SELECTION OF STRATEGY - OUTER-ROI
		/*
		title = "Particle_in_Cell-3D";
		msg = "SEGMENTATION OF THE CELL"
		+"\n \nCheck the outlines and choose a segmentation strategy "
		+"\nto be used (S1 or S2).";
		waitForUser(title, msg);
		do {
			Dialog.create("Particle_in_Cell-3D");
			Dialog.addMessage("SEGMENTATION OF THE CELL");
			Dialog.addMessage("Please choose one option:");
			Dialog.addCheckbox("Strategy S1", false);
			Dialog.addCheckbox("Strategy S2", false);
			Dialog.show();
			S1 = Dialog.getCheckbox();
			S2 = Dialog.getCheckbox();
			S = S1 + S2;
		} while (S != 1);
		*/
		S1 = true;
		if (S1 == true) {
			seg = "S1";
			roiManager("Open", dir_res + "S1-" + cell_name + ".zip");
			roiManager("Deselect");
			roiManager("Save", dir_res + "OUTER-ROI-" + cell_name + ".zip");
			selectWindow("S2-Outlines");
			close();
			selectWindow("S1-Outlines");
		}
		/*
		if (S2 == true) {
			seg = "S2";
			roiManager("Open", dir_res + "S2-" + cell_name + ".zip");
			roiManager("Deselect");
			roiManager("Save", dir_res + "OUTER-ROI-" + cell_name + ".zip");
			selectWindow("S1-Outlines");
			close();
			selectWindow("S2-Outlines");
		}
		*/
		run("Select None");
		saveAs("tiff", dir_res + "CELL_OUTLINES-" + cell_name);
		roiManager("Open", dir_res + "OUTER-ROI-" + cell_name + ".zip");
		roiManager("select", roiManager("count") - 1);
		run("Clear Outside");
		saveAs("tiff", dir_res + "CELL_OUTLINES_Alone-" + cell_name);
		close();
		open(""+ dir_res + "CELL_OUTLINES_Alone-" + cell_name +".tif");
		cell_outlines = File.name;
		roiManager("Reset");
		setPasteMode("Copy");
	
		// 9.7 ENLARGED MEMBRANE REGION
		// This step builds a membrane region of w pixels.
		// 9.7.1 XY-PLANE - INNER-ROI
		if (memb_2D == true) {
			roiManager("Open", dir_res + "OUTER-ROI-" + cell_name + ".zip");
			if (w > 0 ) {
				roiManager("Select", roiManager("count") - 1);
				run("Enlarge...", "enlarge="+-w);
				roiManager("Update");
				run("Draw", "slice");
			}
			roiManager("Deselect");
			roiManager("Save", dir_res + "INNER-ROI-" + cell_name + ".zip");
			roiManager("Reset");
		}
	
		// 9.7.2 FRAME & FINAL CELLULAR OUTLINES
		/* This creates a frame of width = w1 to be used in all images and ROIs
		 * to acount for the membrane region. Segmentation strategy S1 uses a 
		 * frame of 4 pixels; thus w1 is always equal or greater than 4. */
		if (w < 4) { w1 = 4;
		} else { w1 = w;}
		selectWindow(cell_outlines);
		x = getWidth();
		y = getHeight();
		z = nSlices;
		setSlice(slice_max);
		setTool("rectangle");
		makeRectangle(1 + w1, 1 + w1, x - 2 * (1 + w1), y - 2 * (1 + w1));
		run("Clear Outside", "stack");
		run("Draw", "stack");
		run("Select None");
		saveAs("tiff", dir_res + "CELL_OUTLINES-" + cell_name);
		roiManager("Reset");
		run("Close All");
		open(nps_sub);
		makeRectangle(1 + w1, 1 + w1, x - 2 * (1 + w1), y - 2 * (1 + w1));
		run("Clear Outside", "stack");
		run("Select None");
		saveAs("tiff", dir_res + "PARTICLES-" + nps_name);
		setForegroundColor(255, 255, 255);
	  }
	
	  selectWindow("ROI Manager");
	  run("Close");
	  run("Close All");
	  setBatchMode(true);
	  
	  if (R1 || R2 || R3 == true) {
		// 9.7.3 VOLUME AND SURFACE AREA OF THE CELL
		open(cell_mask);
		run("Select All");
		run("Clear", "stack");
		run("8-bit");
		roiManager("Reset");
		roiManager("Open", dir_res + "OUTER-ROI-" + cell_name + ".zip");
		run("Clear Results");
		run("Set Measurements...", "area perimeter stack limit redirect=None decimal=2");
		for(i = 0; i <= (cell_basal - cell_top); i++) {
			roiManager("Select", roiManager("count") - 1);
			run("Fill", "slice");
			makeRectangle(1 + w1, 1 + w1, x - 2 * (1 + w1), y - 2 * (1 + w1));
			run("Clear Outside", "slice");
			setThreshold(1, 255);
			run("Create Selection");
			run("Measure");
			cell_area = cell_area + getResult("Area", i);
			if (i == 0) {
				area_top = getResult("Area", 0);
				area_top_u = (getResult("Area", 0) * xy_scale * xy_scale) / 1E6;
			} else {
				area_diff = abs(getResult("Area", i) - getResult("Area", i - 1));
				area_ratio = area_diff / getResult("Perim.", i);
				if (area_ratio < 0.90) {
					area = area + (getResult("Perim.", i));
					area_u = area_u + ((getResult("Perim.", i) * z_scale * xy_scale) / 1E6);
				}
				if (area_ratio >= 0.90 && area_ratio < 1.11) {
					area = area + ((getResult("Perim.", i) + area_diff));
					area_u = area_u + ((getResult("Perim.", i) * z_scale * xy_scale) / 1E6) + 
					((area_diff * xy_scale * xy_scale) / 1E6);
				}
				if (area_ratio >= 1.11) {
					area = area + area_diff;
					area_u = area_u + ((area_diff * xy_scale * xy_scale) / 1E6);
				}
				if (i == (cell_basal - cell_top)) {
					area_basal = getResult("Area", i);
					area_basal_u = (area_basal * xy_scale * xy_scale) / 1E6;
				}
			}
		}
		area_apical = area_top + area;
		area_apical_u = area_top_u + area_u;
		cell_area_u = (area_unit * cell_area)/1E9;
		run("Select None");
		saveAs("tiff", dir_res + "CELL_MASK_VOLUME-" + cell_name);
		roiManager("Reset");
		run("Clear Results");
	
		// 9.7.4 INTRACELLULAR VOLUME - MEMBRANE REGION - XY-PLANE
		if (memb_2D == true) {
			roiManager("Open", dir_res + "INNER-ROI-" + cell_name + ".zip");
			run("Set Measurements...", "area limit redirect=None decimal=0");
			for (i = 0; i < w_slices; i++) {
				roiManager("Select", roiManager("count") - 1);
				run("Select All");
				run("Clear", "slice");
			}
			for (i = w_slices; i <= (cell_basal - cell_top); i++) {
				roiManager("Select", roiManager("count") - 1);
				run("Select All");
				run("Clear", "slice");
				roiManager("Select", roiManager("count") - 1);
				run("Fill", "slice");
				setThreshold(1, 255);
				run("Measure");
				cell_area_in = cell_area_in + getResult("Area", i - w_slices);
			}
			run("Select None");
			saveAs("tiff", dir_res + "CELL_MASK_INTRACELL-" + cell_name);
		}
		cell_area_in_u = (area_unit * cell_area_in)/1E9;
		run("Clear Results");
		run("Close All");
	
	//NEW	
		// 9.7.5 MEMBRANE REGION - MASK - XY-PLANE
		open(""+ dir_res + "CELL_MASK_VOLUME-" + cell_name +".tif");
		if (memb_2D == true) {
			for (i = w_slices; i <= (cell_basal - cell_top); i++) {
				roiManager("Select", roiManager("count") - 1);
				run("Clear", "slice");
			}
			run("Select None");
			saveAs("tiff", dir_res + "CELL_MASK_MEMBRANE-" + cell_name);
		}
		run("Clear Results");
		roiManager("Reset");
		run("Close All");
	
	//END NEW
	  }
	// 10 EVALUATION OF PARTICLES
	// 10.1 SEGMENTATION OF PARTICLES WITH PLUGIN '3D OBJECT COUNTER'
	  open(nps_seg);
	  run("Select None");
	  if (imj_fiji == true) {
		run("3D OC Options", "volume surface nb_of_obj._voxels integrated_density "
		+"mean_gray_value maximum_gray_value centre_of_mass bounding_box "
		+"close_original_images_while_processing_(saves_memory) "
		+"redirect_to=none");
		if (exclude == true) {
			run("3D Objects Counter", "threshold="+thd_lower+" slice="+slice_max+" "
			+"min.="+obj_min_area+" max.="+obj_max_area+" exclude_objects_on_edges "
			+"objects statistics");
		}
		if (exclude == false) {
			run("3D Objects Counter", "threshold="+thd_lower+" slice="+slice_max+" "
			+"min.="+obj_min_area+" max.="+obj_max_area+" objects statistics");
		}
	  }
	  if (imj == true) {
		run("3D Manager Options", "volume nb_of_obj._voxels integrated_density "
		+"mean_gray_value maximum_gray_value centre_of_mass bounding_box "
		+"close_original_images_while_processing_(saves_memory) "
		+"redirect_to=none");
		if (exclude == true) {
			run("Object Counter3D", "threshold="+thd_lower+" slice="+slice_max+" "
			+"min.="+obj_min_area+" max.="+obj_max_area+" exclude_objects_on_edges "
			+"objects statistics");
		}
		if (exclude == false) {
			run("Object Counter3D", "threshold="+thd_lower+" slice="+slice_max+" "
			+"min.="+obj_min_vol+" max.="+obj_max_vol+" objects statistics");
		}
	  }
	  run("8-bit");
	  for (i = 0; i < nSlices; i++) {
		setSlice(i + 1);
		changeValues(1,254,255);
	  }
	  run("Grays");
	  run("RGB Color");
	  saveAs("tiff", dir_res + "PARTICLES_SEGMENTED-" + nps_name);
	  run("Close All");
	
	// 10.2 POSITION OF OBJECTS AND BOUNDING BOXES
	  obj_nr = nResults;
	  bx = newArray(obj_nr); // Bounding boxes
	  by = newArray(obj_nr);
	  bz = newArray(obj_nr);
	  bxx = newArray(obj_nr);
	  byy = newArray(obj_nr);
	  bzz = newArray(obj_nr);
	  obj_x = newArray(obj_nr); // X
	  obj_y = newArray(obj_nr); // Y
	  obj_z = newArray(obj_nr); // Z
	  obj_vol = newArray(obj_nr); // Volume
	  if (MultCells) {
	  	CellNum = newArray(obj_nr); // Cell from wich they come	
	  }
	  
	  for (i = 0; i < obj_nr; i++) {
		bx[i] = getResult("BX", i);
		by[i] = getResult("BY", i);
		bz[i] = getResult("BZ", i);
		bxx[i] = getResult("B-width", i);
		byy[i] = getResult("B-height", i);
		bzz[i] = getResult("B-depth", i);
		obj_x[i] = getResult("XM", i);
		obj_y[i] = getResult("YM", i);
		obj_z[i] = round(getResult("ZM", i));
		obj_vol[i] = getResult("Nb of obj. voxels", i);
		if (MultCells) {
			for (j = 0; j < Roix.length; j++){
				if (obj_x[i] > Roix[j] && obj_x[i] < Roix[j] + RoiWidth[j] && obj_y[i] > Roiy[j] && obj_y[i] < Roiy[j] + RoiHeight[j]) {
					CellNum[i] = j + 1;
					break; 
				}
				else {
					CellNum[i] = 0;
				}
			}
		}		
		if (obj_z[i] == 0) {
			obj_z[i] = 1; }
		if (obj_z[i] == nps_bottom + 1) {
			obj_z[i] = nps_bottom;}
	  }
	  selectWindow("Results");
	  saveAs("Results", dir_res + "TABLE_3DOC-" + nps_name + ".xls");
	  run("Close");
	
	// 10.3 INTENSITY AND NUMBER OF PARTICLES IN EACH OBJECT
	  open(nps_sub);
	  obj_id = newArray(obj_nr);
	  obj_nr_nps = newArray(obj_nr);
	  run("Set Measurements...", "  integrated redirect=None decimal=0");
	  for (i = 0; i < obj_nr; i++) {
		run("Select None");
		makeRectangle(bx[i], by[i], bxx[i], byy[i]);
		run("Measure");
		obj_id[i] = obj_id[i] + getResult("RawIntDen", 0);
		obj_nr_nps[i] = round(obj_id[i] / mean_id);
		run("Clear Results");
	  }
	  run("Select None");
	  if (R1 || R2 || R3 == true) {
	
		// 10.4 LOCATION AND NUMBER OF PARTICLES INTERACTING WITH THE CELL
		obj_in = newArray(obj_nr);
		obj_memb = newArray(obj_nr);
		obj_memb_basal = newArray(obj_nr);
		obj_out = newArray(obj_nr);
		run("Duplicate...", "title=Locator");
		selectWindow(nps_sub);
		close();
		selectWindow("Locator");
		run("8-bit");
		run("Select All");
		run("Clear", "slice");
		if (memb_2D == true) {
			roiManager("Open", dir_res + "INNER-ROI-" + cell_name + ".zip");
	}
		roiManager("Open", dir_res + "OUTER-ROI-" + cell_name + ".zip");
		run("Set Measurements...", "area redirect=None decimal=0");
		run("Clear Results");
		for (i = 0; i < obj_nr; i++) {
			locator = 0;
			selectWindow("Locator");
			run("Select All");
			run("Clear", "slice");
			run("Select None");
			setPixel(obj_x[i], obj_y[i], 255);
	
			// 10.4.1 INTRACELLULAR PARTICLES
			if (memb_2D == true) {
				if (obj_z[i] >= cell_top + w_slices && obj_z[i] <= (nps_bottom)) {
					roiManager("Select", roiManager("count") - 1);
					run("Analyze Particles...", "size=0-1 pixel"
					+" circularity=0.00-1.00 show=Nothing clear");
					locator = nResults;
					if (locator == 1) {
						obj_in[i] = obj_nr_nps[i];
						obj_memb[i] = 0;
					}
					if (locator == 0) {
						obj_in[i] = 0;
					}
				}
			}
			if (memb_3D == true) {
				if (obj_z[i] >= cell_top + w_slices && obj_z[i] <= (cell_basal_end)) {
					roiManager("Select", obj_z[i] - cell_top);
					run("Clear", "slice");
					roiManager("Select", cell_basal_end - cell_top + 1 + obj_z[i] - cell_top);
					run("Clear Outside", "slice");
					run("Select All");
					run("Analyze Particles...", "size=0-1 pixel"
					+" circularity=0.00-1.00 show=Nothing clear");
					locator = nResults;
					if (locator == 1) {
						obj_in[i] = obj_nr_nps[i];
						obj_memb[i] = 0;
						obj_memb_basal[i] = 0;
					}
					if (locator == 0) {
						obj_in[i] = 0;
					}
				}
			}
	
			// 10.4.2 MEMBRANE ASSOCIATED PARTICLES 
			if (locator == 0) {
				
	if (memb_2D == true) {
					roiManager("Select", nps_bottom - cell_top + 1 + obj_z[i] - cell_top);
					run("Analyze Particles...", "size=0-1 pixel"
					+" circularity=0.00-1.00 show=Nothing clear");
					locator = nResults;
					if (locator == 1) {
						obj_memb[i] = obj_nr_nps[i];
						obj_in[i] = 0;
					}
					if (locator == 0) {
						obj_memb[i] = 0;
					}
				}
				if (memb_3D == true) {
					setPixel(obj_x[i], obj_y[i], 255);
					roiManager("Select", cell_basal_end - cell_top + 1 + obj_z[i] - cell_top);
					run("Analyze Particles...", "size=0-1 pixel"
					+" circularity=0.00-1.00 show=Nothing clear");
					locator = nResults;
					if (locator == 1) {
						if (obj_z[i] >= (cell_basal-round(w_slices/2)+1) && obj_z[i]<=(cell_basal_end)) {
							obj_memb_basal[i] = obj_nr_nps[i];
							obj_memb[i] = 0;
							obj_in[i] = 0;
						} else {
							obj_memb[i] = obj_nr_nps[i];
							obj_memb_basal[i] = 0;
							obj_in[i] = 0;
							if (w == 0) {
								obj_in[i] = obj_nr_nps[i];
								obj_memb[i] = 0;
								obj_memb_basal[i] = 0;
							}
						}
					}
					if (locator == 0) {
						obj_memb[i] = 0;
						obj_memb_basal[i] = 0;
					}
				}
			}
			// 10.4.3 EXTRACELLULAR PARTICLES 
			if (locator == 0) {
				obj_out[i] = obj_nr_nps[i];
	}
			if ((obj_z[i] >= 1 && obj_z[i] < cell_top) || (obj_z[i] > cell_basal_end)) {
				obj_out[i] = obj_nr_nps[i];
				obj_memb[i] = 0;
				obj_memb_basal[i] = 0;
			}
		}
		run("Clear Results");
		selectWindow("Locator");
		close();
		roiManager("Reset");
	  }
	
	// 10.5 RESULTS TABLES
	  obj_c = 0;
	  obj_c_in = 0;
	  obj_c_memb = 0;
	  obj_c_memb_basal = 0;
	  total_id_in = 0;
	  total_id_memb = 0;
	  total_id_memb_basal = 0;
	  if (R1 || R2 || R3 == true) {
		for (i = 0; i < obj_nr; i++) {
			if (obj_in[i] >= 1 || obj_memb[i] >= 1 || obj_memb_basal[i] >= 1) {
				setResult("3DOC Label", obj_c, i + 1);
				setResult("Object", obj_c, obj_c + 1);
				setResult("Volume (vx)", obj_c, obj_vol[i]);
				setResult("IntDens th=0 (px value)", obj_c, obj_id[i]);
				setResult("Particles Intracellular", obj_c, obj_in[i]);
				if (memb_2D == true) {
					setResult("Particles Membrane", obj_c, obj_memb[i]);
	}
				if (memb_3D == true) {
					setResult("Particles Apical Memb.", obj_c, obj_memb[i]);
					setResult("Particles Basal Memb.", obj_c, obj_memb_basal[i]);
				}
				setResult("Particles Total", obj_c, obj_nr_nps[i]);
				setResult("X (px)", obj_c, obj_x[i]);
				setResult("Y (px)", obj_c, obj_y[i]);
				setResult("Z (slice)", obj_c, obj_z[i]);
				if (MultCells) {
					setResult("Cell", obj_c, CellNum[i]);	
				}
				
				if (obj_in[i] >= 1 ) {
					obj_c_in = obj_c_in + 1;
					total_id_in = total_id_in + obj_id[i];
				}
				if (obj_memb[i] >= 1 ) {
					obj_c_memb = obj_c_memb + 1;
					total_id_memb = total_id_memb + obj_id[i];
				}
				if (memb_3D == true) {
					if (obj_memb_basal[i] >= 1 ) {
						obj_c_memb_basal = obj_c_memb_basal + 1;
						total_id_memb_basal = total_id_memb_basal + obj_id[i];
					}
				}
				obj_c = obj_c + 1;
			}
		}
		for (i = 0; i < obj_c; i++) {
		nps_in = nps_in + getResult("Particles Intracellular", i);
			if (memb_2D == true) {
				nps_memb = nps_memb + getResult("Particles Membrane", i);}
			if (memb_3D == true) {
				nps_memb = nps_memb + getResult("Particles Apical Memb.", i);
				nps_memb_basal = nps_memb_basal + getResult("Particles Basal Memb.", i);
			}
		}
		if (memb_2D == true) {
			total_id = total_id_in + total_id_memb;
			nps_total = nps_in + nps_memb;
		}
		if (memb_3D == true) {
			total_id = total_id_in + total_id_memb + total_id_memb_basal;
			nps_total = nps_in + nps_memb + nps_memb_basal;
		}
	  }
	  if (R4 || R5 == true) {
		for (i = 0; i < obj_nr; i++) {
			setResult("3DOC Label", obj_c, i + 1);
			setResult("Object", obj_c, obj_c + 1);
			setResult("Volume (vx)", obj_c, obj_vol[i]);
			setResult("IntDens th=0 (px value)", obj_c, obj_id[i]);
			if (R5 == true) {
			setResult("Number of Particles", obj_c, obj_nr_nps[i]);}
			setResult("X (px)", obj_c, obj_x[i]);
			setResult("Y (px)", obj_c, obj_y[i]);
			setResult("Z (slice)", obj_c, obj_z[i]);
			total_id = total_id + obj_id[i];
			obj_c = obj_c + 1;
		}
	  }
	  if (R5 == true) {
	  	for (i = 0; i < obj_c; i++) {
			nps_total = nps_total + getResult("Number of Particles", i);
		}
	  }
	  setOption("ShowRowNumbers", false);
	  updateResults;
	  saveAs("Results", dir_res + "_RESULTS-" + nps_name + ".xls");
	  saveAs("Results", dir_res + "_RESULTS-" + nps_name + ".txt");
	  if (R1 || R2 || R3 == true) {
	
		// 10.6 COLOR CODING OF PARTICLES
		open(nps_seg);
		for (i = 0; i < obj_nr; i++) {
			// 10.6.1 INTRACELLULAR
			if (obj_in[i] >= 1) {
				makeRectangle(bx[i], by[i], bxx[i], byy[i]);
					if (color_nps_in == "Red") {
						changeValues(0xffffff,0xffffff,0xff0000);}
					if (color_nps_in == "Green") {
						changeValues(0xffffff,0xffffff,0x00ff00);}
					if (color_nps_in == "Magenta") {
						changeValues(0xffffff,0xffffff,0xff00ff);}
			}
	
			// 10.6.2 MEMBRANE ASSOCIATED
			// 10.6.2.1 MEMBRANE ASSOCIATED - APICAL
			if (obj_memb[i] >= 1) {
				makeRectangle(bx[i], by[i], bxx[i], byy[i]);
					if (color_nps_memb == "Yellow") {
						changeValues(0xffffff,0xffffff,0xffff00);}
					if (color_nps_memb == "Cyan") {
						changeValues(0xffffff,0xffffff,0x00ffff);}
			}
	
	
			// 10.6.2.2 MEMBRANE ASSOCIATED - BASAL
			if (memb_3D == true) {
				if (obj_memb_basal[i] >= 1) {
					makeRectangle(bx[i], by[i], bxx[i], byy[i]);
					for (k = bz[i]; k < (bz[i] + bzz[i]); k++) {
						setSlice(k);
						if (color_nps_memb_basal == "Orange") {
							changeValues(0xffffff,0xffffff,0xffb900);}
						if (color_nps_memb_basal == "Yellow") {
							changeValues(0xffffff,0xffffff,0xffff00);}
						if (color_nps_memb_basal == "Cyan") {
							changeValues(0xffffff,0xffffff,0x00ffff);}
						if (color_nps_memb_basal == "Red") {
							changeValues(0xffffff,0xffffff,0xff0000);}
						if (color_nps_memb_basal == "Green") {
							changeValues(0xffffff,0xffffff,0x00ff00);}
						if (color_nps_memb_basal == "Magenta") {
							changeValues(0xffffff,0xffffff,0xff00ff);}
					}
				}
			}
	
			// 10.6.3 EXTRACELLULAR
			if (obj_out[i] >= 1) {
				makeRectangle(bx[i], by[i], bxx[i], byy[i]);
				for (k = bz[i]; k < (bz[i] + bzz[i]); k++) {
					setSlice(k);
					if (color_nps_out == "Do Not Show") {
						changeValues(0xffffff,0xffffff,0x000000);}
				}
			}
	
			// 10.6.4 CLEARING OBJECTS FORMED BY "ZERO" PARTICLES
			if (obj_nr_nps[i] == 0) {
				makeRectangle(bx[i], by[i], bxx[i], byy[i]);
				for (k = bz[i]; k < (bz[i] + bzz[i]); k++) {
					setSlice(k);
					changeValues(0xffffff,0xffffff,0x000000);
				}
			}
		}
	
		// 10.7 COLORED VISUALIZATION OF PARTICLE UPTAKE IN 3D
		
		// 10.7.1 Z-PROJECTION
		selectWindow(nps_seg);
		run("Select None");
		run("Gaussian Blur...", "sigma=1.0 stack");
		wait(1000);
		saveAs("tiff", dir_res + "PARTICLES_SMOOTHED-" + nps_name);
		roiManager("Open", dir_res + "OUTER-ROI-" + cell_name + ".zip");
	 	roiManager("Select", roiManager("count") - 1);
	 	run("Clear Outside");
	 	saveAs("tiff", dir_res + "PARTICLES_SMOOTHED_Alone-" + nps_name);
		close();
		open(""+ dir_res + "PARTICLES_SMOOTHED_Alone-" + nps_name +".tif");
		nps_smo = File.name;
		if (memb_2D == true) {
			open(""+ dir_res + "CELL_MASK_MEMBRANE-" + cell_name +".tif");}
		rename("Membrane");
		setAutoThreshold("IJ_IsoData dark");
		//doWand(xSeed[0], ySeed[0]);
		
		// 10.7.2 UPTAKE - CELL MASK
		selectWindow("Membrane");
		run("8-bit");
		run(color_cell);
		run("RGB Color");
		imageCalculator("Transparent-zero create", ""
		+"Membrane",""+nps_smo+"");
		saveAs("tiff", dir_res + "_UPTAKE_MASK-" + cell_name);
		
		// 10.7.3 UPTAKE
		open(cell_outlines);
		imageCalculator("Transparent-zero create stack", ""
		+""+cell_outlines+"",""+nps_smo+"");
		saveAs("tiff", dir_res + "_UPTAKE-" + cell_name);
		close();
		open(""+ dir_res + "_UPTAKE-" + cell_name +".tif");
		nps_uptake = File.name;
		run("Close All");
	  }
	
		// 10.8 COLOR PARTICLES - CALIBRATION & ONLY PARTICLES 
	  if (R4 || R5 == true) {  
		open(nps_seg);
		run("8-bit");
		run(""+color_nps_in+"");
		run("RGB Color");
		run("Gaussian Blur...", "sigma=1.0 stack");
		setMinAndMax(0, 255);
		saveAs("tiff", dir_res + "PARTICLES_SMOOTHED-" + nps_name);
		close();
		open(""+ dir_res + "PARTICLES_SMOOTHED-" + nps_name +".tif");
		nps_smo = File.name;
		run("Close All");
		selectWindow("Results");
		run("Close");
	  }
	
	// 11 REPORTS
	  setBatchMode(false);
	  print("");
	  print("************************************************");
	  print("Particle_in_Cell-3D - ImageJ Macro");
	  print("************************************************");
	  
	print("");
	  print("-------------------------------------------------------------------------------");
	  if (R1 == true) {print("  REPORT OF RESULTS - 1. QUALITATIVE");}
	  if (R2 == true) {print("  REPORT OF RESULTS - 2. SEMI-QUANTITATIVE");}
	  if (R3 == true) {print("  REPORT OF RESULTS - 3. QUANTITATIVE");}
	  if (R4 == true) {print("  REPORT OF RESULTS - 4. CALIBRATION");}
	  if (R5 == true) {print("  REPORT OF RESULTS - 5. ONLY PARTICLES");}
	  print("-------------------------------------------------------------------------------");
	  print("    Experiment: 	" +exp_name);
	  getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
	  t_hour = hour;
	  t_minute = minute;
	  if (hour < 10) {t_hour = "0" +hour; }
	  if (minute < 10) {t_minute = "0" +minute; }
	  print("    Date & Time: 	" +dayOfMonth+ "/" +month+1+ "/" +year+ ", " +t_hour+ ":"
	  +t_minute);
	  print("    Original file-PARTICLES: 	" +nps_path);
	  print("        Substack: 	"+"(" +nps_top_p+ "-" +nps_bottom_p+ ")");
	  if (R1 || R2 || R3 == true) {
		print("    Original file-CELL: 	" +cell_path);
		print("        Substack: 	"+"(" +nps_top_p+ "-" +nps_bottom_p+ ")");
		print("        Selected Z-position of the cell: 	"+"(" +cell_top_p+ "-" +cell_basal_p+ ")");
	  }
	  	print("    Directory with Results: 	"+ dir_res +"");
	  print("");
	  print("-------------------------------------------------------------------------------");
	  print("ANALYSIS PARAMETERS");
	  print("-------------------------------------------------------------------------------");
	  if (R1 == false) {
		print("    XY-scale (nm/px): 	" +xy_scale);
		print("    Voxel size (nm^2): 	" +d2s(area_unit, 0));
	  }
	  print("    Subtracted background (px value): 	" +bg);
	  print("    Threshold - Particles (px value): 	" +thd_lower);
	  if (R3 || R5 == true) {
		print("    Mean intensity of a single particle (px value): 	" +mean_id);
	  }
	  print("    Include particles larger than (voxels): 	" +obj_min_area);
	  print("        and smaller than (voxels): 	" +obj_max_area);
	  if (exclude == false) {
	  	print("    Exclude objects on edges: 	yes");}
	
	  if (exclude == false) {
	  	print("    Exclude objects on edges: 	no");}
	  if (R1 || R2 || R3 == true) {
	  	print("");
	  	print("    Threshold - Cell (px value): 	" +thd_c_lower);
	  	if (M1 == true) {
	  		print("    Mask-1 (Automatic threshold)");
	  	}
	  	if (M2 == true) {
	  		print("    Mask-2");
	  	}
		print("    Segmentation strategy: 	" +seg);
		if (memb_2D == true) {
			print("    Membrane Region XY-Plane:");
			print("    Width (px, XY-scale): 	" +w);
			print("    Width (�m): 	" +d2s(((w*xy_scale)/1000), 3));
		}
		if (memb_3D == true) {
			print("    Membrane Region-3D:");
			print("    Width (px, XY-scale): 	" +w);
			print("    Width (�m): 	" +d2s(((w*xy_scale)/1000), 3));
			print("    Off-set (px, XY-scale): 	" +w_offset);
		}
		
		if (PreciseQuant) {
			print("    Precise Quantification Parameters 	");
			print("    Flourescence Ratio: 	" + FluorescenseRatio);
			print("    Calibration Factor: 	" + CalibrationFactor);
		}
	  }
	  print("");
	  print("-------------------------------------------------------------------------------");
	  if (R1 == false) {
		print("RESULTS");
		print("-------------------------------------------------------------------------------");
		if (R2 || R3 == true) {
			print("    CELL:");
			print("        Area of basal membrane (px): 	" +d2s(area_top, 1));
			print("        Area of basal membrane (�m^2): 	" +d2s(area_top_u, 1));
			print("        Area of apical membrane (px): 	" +d2s(area_apical, 1));
			print("        Area of apical membrane (�m^2): 	" +d2s(area_apical_u, 1));
			print("    NUMBER OF OBJECTS:");
			print("        Intracellular: 	" +obj_c_in);
			if (memb_2D == true) {
				print("        Membrane region: 	" +obj_c_memb);}
			print("        Total: 	" +obj_c);
		}
		if (R3 == true) {
			print("    NUMBER OF PARTICLES:");
			print("        Intracellular: 	" +nps_in);
			if (memb_2D == true) {
				print("        Membrane region: 	" +nps_memb);}
			print("        Total: 	" +nps_total);
			print("    DENSITY OF PARTICLES:");
			print("        Intracellular (particles/1000*nm^3): 	" +d2s(1000 * nps_in / cell_area_in_u, 1));
			print("        Whole cell (particles/1000*nm^3): 	" +d2s(1000 * nps_total / cell_area_u, 1));
			print("        Apical membrane (particles/100*nm^2): 	" +d2s(100 * nps_memb / area_apical_u, 1));
		}
		if (R2 || R3 == true) {
			print("    PIXEL INTENSITY OF PARTICLES:");
			print("        Intracellular: 	" +total_id_in);
			if (memb_2D == true) {
				print("        Membrane region: 	" +total_id_memb);}
			print("        Total: 	" +total_id);
		}
		if (R4 || R5 == true) {
			print("    Number of Objects: 	" +obj_c+"");}
		if (R5 == true) {
			print("    Pixel intensity (total): 	" +total_id);
			print("    Number of particles: 	" +nps_total);
			print("    Area analyzed (voxels): 	" +cell_area);
			print("    area analyzed (nm^3): 	" +d2s(cell_area_u, 1));
			print("    Number density (particles/1000*nm^3): 	" +d2s(1000 * nps_total / cell_area_u, 1));
		}
	
	
		
		print("");
		print("-------------------------------------------------------------------------------");
		print("TABLE WITH DETAILED RESULTS AT:");
		print("-------------------------------------------------------------------------------");
		print("    File: 	_RESULTS-"+ nps_name + ".xls");
		print("    Directory: 	"+ dir_res +"");
		print("");
		print("-------------------------------------------------------------------------------");
	  }
	  print("	***End***");
	
	  // 11.1 SAVING FILES AND VISUALIZATION OF RESULTS
	  // Saving the Report of Results
	  selectWindow("Log");
	  saveAs("text", dir_res + "_ANALYSIS_REPORT-" + nps_name + ".txt");
	  selectWindow("Log");
	  run("Close");
	  run("Text File... ", "open=[" + dir_res + "_ANALYSIS_REPORT-" + nps_name + ".txt]");
	  if (R1 || R2 || R3 == true) {
	  	open(nps_uptake);
		setSlice(slice_max);
	  }
	  if (R4 || R5 == true) {
	  	open(nps_smo);
		setSlice(slice_max);
		open(nps_seg);
		setSlice(slice_max);
		run("Tile");
	  }
	  if (MultIm) {
	  	Dialog.create("Particles in Cell");
	  	Dialog.addCheckbox("More Images?", false);
	  	Dialog.show();
	  	MultIm = Dialog.getCheckbox();
	  }
} while(MultIm);

exit("Particle_in_Cell-3D\n \nProcessing is finished!")