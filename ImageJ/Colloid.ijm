// choose project folder

Dir = getDirectory("Choose analysis folder:");

//Image Directory
Dir_im = ""+Dir+"/Images/";


//Create folder for 'Completed Images'
	// Dir_good
	// dump roi and images here
Dir_good = ""+Dir+"/CompletedImages/";
File.makeDirectory(Dir_good);

//Create folder for 'NeedToLookAt'
	//Dir_nogood
Dir_nogood = ""+Dir+"/NeedToLookAt/";
File.makeDirectory(Dir_nogood);

//Create folder for .csv
	//Dir_csv
Dir_csv = ""+Dir+"/CSV/";
File.makeDirectory(Dir_csv);
//
//
//-----------------------------------------------
run("Set Measurements...", "area mean min integrated display redirect=None decimal=3");
//------------------------------------------------
//
//
//Functions --------------------------------------
//
	//---------------------------
	//----SelectCheck------------
	//---------------------------
	// Can save the program from crashing if you're supposed
	//   to have an area selected. If there is no area selected, 
	//   the program may exit. This runs a check to make sure 
	//	 SOMETHING is selected before moving forward and potentially
	//	 crashing the program. 
	
	function selectCheck() {
		checkVal = getValue("selection.size");
		if (checkVal == 0){
			waitForUser("You forgot to check select!");
		}
	}

	//---------------------------
	//----SelectAllSelectNone----
	//---------------------------
	// Quick way to remove the roi selections, which can
	//   stuff up masking and relevant.
	
	function noSelect() {
	 	roiManager("Show All");
		roiManager("Show None");
	}

	//---------------------------

//------------------------------------------------
//
// Roi manager
// [0] = Colloids
// [1] = Colloids Threshold
//
//------------------------------------------------
//
nextImage = "Yes";
while(nextImage =="Yes"){
	//open image
		//get file name
	images = getFileList(Dir_im);
	len = images.length;
	if (len >1 ){
		randnum = Math.floor(random()*len);
	} else {
		randnum =1;
	}
	//randnum = 0; //uncomment and comment above to un-randomise
	workingImage = images[randnum];
	open(Dir_im+workingImage);


		// First opportunity to scrap it. 
	Dialog.create("Is this image Good?");
	Dialog.addChoice("Are you happy to continue? Selecting No will move this to the 'no good' folder.", newArray("Yes", "No"), "Yes");
	Dialog.show();
	moveImage = Dialog.getChoice();

		// if no
			// move image to Dir_nogood
		if(moveImage=="No") {
			newDir = Dir_nogood+workingImage;
			File.rename(Dir_im+workingImage,newDir);
			close("*");
		}// end no
		
	// else continue
	else {
		//Draw around your colloids
		setTool("polygon");
		waitForUser("Outline your colloid area\n \nBe detailed!\n \nDo not click OK until outline complete!");
		selectCheck();
		roiManager("Add");
		roiManager("Select", 0);
		roiManager("Rename", "Colloids");
		noSelect();
		
			//open Color Threshold and prompt the user on how to threshold.
		run("Color Threshold...");
		waitForUser("Threshold the whiteness using Brightness.\n \nAre there breaks in your outlines?\n \nClick Select:\n \nDo not click OK until you have selected!");
		selectCheck();
		selectCheck();
		roiManager("Add");
		roiManager("Select", 1);
		roiManager("Rename", "Colloid Threshold");
	
			//make Colloid mask
		roiManager("Select", 1);
		run("Create Mask");
		roiManager("Select", 0);
		run("Make Inverse");
		setBackgroundColor(0, 0, 0);
		run("Clear", "slice");
			
			// Get min size
		selectWindow(workingImage);
		setTool("oval");
		noSelect();
		waitForUser("Outline the INSIDE of the smallest colloid\n \nDo not click OK until you have done so!");	
		run("Measure");
		selectWindow("Results");
		small = getResult("Area", 0);
		Table.deleteRows(0, 0);
	
			// Get max size
		selectWindow(workingImage);
		setTool("oval");
		noSelect();
		waitForUser("Outline the OUTSIDE of the largest colloid\n \nDo not click OK until you have done so!");
		run("Measure");
		selectWindow("Results");
		large = getResult("Area", 0);
		Table.deleteRows(0, 0);
	
			// Make mask containing only necessary colloids
		selectWindow("Mask");
		noSelect();
		run("Analyze Particles...", "size=small-large show=Masks display summarize");
		close("Results");
	
			// Opportunity to remove unwanted spots
		selectWindow("Mask of Mask");
		run("Invert");
		setTool("oval");
		selectWindow(workingImage);
		run("Out [-]");run("Out [-]");
		selectWindow("Mask of Mask");
		waitForUser("Do not press continue until you have checked this Mask with the\noriginal image and removed unwanted spots.");
	
			// Second opportunity to scrap it. 
		Dialog.create("Is this image Good?");
		Dialog.addChoice("Are you happy with the Analysis??\n \nSelecting No will move this to the 'no good' folder.", newArray("Yes", "No"), "Yes");
		Dialog.show();
		moveImage = Dialog.getChoice();
	
			// if no
				// move image to Dir_nogood
			if(moveImage=="No") {
				newDir = Dir_nogood+workingImage;
				File.rename(Dir_im+workingImage,newDir);
			}// end no
		
		//else contonue
	
		noSelect();
	
		
			//Check Distance Map
		close("Mask");
		selectWindow("Mask of Mask");
		run("Invert");
		run("Distance Map");
		run("Threshold...");
		setThreshold(0, 20);
		
			//Ask: Correct orientation? 
		Dialog.create("Is this the correct orientation for the distance map?");
			Dialog.addChoice("", newArray("Yes", "No")) 
			Dialog.show();
			choice2 = Dialog.getChoice();
		if (choice2 == "No") {
			run("Undo");
			resetThreshold();
			selectWindow("Mask of Mask");
			run("Invert");
			run("Distance Map");
			run("Threshold...");
			setThreshold(0, 20);
		} //end if orientation is correct
	
			//Average epithelium
		selectWindow(workingImage);
		setTool(4);
		waitForUser("Measure Epithlium at least 10x");	
		selectCheck();
		arr = newArray(getValue("results.count")); 
		for (i =0; i<getValue("results.count"); i++){
			arr[i] = getResult("Length", i);
	 	}
		Array.getStatistics(arr,nil1, nil2, mean);
		close("Results");
	
			// Set Threshold
		selectWindow("Mask of Mask");
		//noSelect();
		run("Threshold...");
		noSelect();
		setThreshold(1, mean);
		run("Create Selection");
		roiManager("Add");
		roiManager("Select", 2);
		roiManager("Rename", "Epithelium");
	
			// Run Colour Deconv
		selectWindow(workingImage);
		run("Colour Deconvolution", "vectors=[H DAB] hide");
		colour1 = ""+workingImage+"-(Colour_1)";
		colour2 = ""+workingImage+"-(Colour_2)";
		colour3 = ""+workingImage+"-(Colour_3)";
		close(colour3);close(colour1);
		selectWindow(colour2);
		run("Invert");
	
		noSelect();
		roiManager("Select", 0);
		run("Measure");
		roiManager("Select", 1);
		run("Measure");
		roiManager("Select", 2);
		run("Measure");
	
			// Final opportunity to scrap it. 
		Dialog.create("Is this image Good?");
		Dialog.addChoice("Are you happy with the Analysis??\n \nSelecting No will move this to the 'no good' folder.", newArray("Yes", "No"), "Yes");
		Dialog.show();
		
		//Are you happy? 
				Dialog.create("Is this image Good?");
				Dialog.addChoice("Are you happy with the thresholding??\n \nSelecting No will move this to the 'no good' folder.", newArray("Yes", "No"), "Yes");
				Dialog.show();
				moveImage = Dialog.getChoice();
	
			// if no
				// move image to Dir_nogood
				if(moveImage=="No") {
					newDir = Dir_nogood+workingImage;
					File.rename(Dir_im+workingImage,newDir);
				}
			//if yes	
				//Math and saving everything
				if(moveImage == "Yes"){
					
					newDir = Dir_good+workingImage;
					File.rename(Dir_im+workingImage,newDir);
					selectWindow("Results");
					saveAs("Results", ""+Dir_csv+workingImage+"_StainResults.csv");
					noSelect();
					roiManager("Save", ""+Dir_good + workingImage +"_ROIs.zip");
				}
				close("*");
				close("Results");
				roiManager("deselect");
				roiManager("delete");
	}			
	Dialog.create("Do you want to open another image?");
	Dialog.addChoice("If you select yes, it will select the next image in the list.\n \nIf you select No, this program will end.", newArray("Yes", "No"), "Yes");
	Dialog.show();

	nextImage = Dialog.getChoice();
}//End While loop
//