---> Workflow for the code

Global selection:

		Enter component thickness in cm (variable b in diffusion-modified and d in creep-modified)
		Enter temperature in degree Celsius (variable Temp in diffusion-modified and temp in creep-modified)
		Enter stress in MPa (variable stress in creep-modified)
		Enter end time for simulation (variable final_time in diffusion-modified)

1. User chooses materials:

	- Options are (can be drop-down menu)
		a. 740H
		b. 282
		c. 625 
		

		
---> For corrosion module:

2. User chooses environment

	- Options are (can be drop-down menu)
		a. sCO2
		b. KCl-MgCl2
		
3. Based on material selections:

	- if mat=740H
		read input file CorrProp_740H.txt
	  elseif mat=282
		read input file CorrProp_282.txt
	  elseif mat=625
		read input file CorrProp_625.txt
	  end
	  
4. Based on environment

	- if env=sCO2
		set corr_var=1
	  else
		set corr_var=2
	  end
	  
5. Output is :

	if corr_var=1

		- time to 5 wt.% concentration (if variable surfaceconcentration<0.05 generate message "corrosion lifetime criteria met")
	else
		- depth of attack (variable attack)
	
	
---> For creep module:

6. Based on material selections:

	- if mat=740H
		read input file CreepProp_740H.txt
	  elseif mat=282
		read input file CreepProp_282.txt
	  elseif mat=625
		read input file CreepProp_625.txt
	  end
	  
7. Output is (I have marked these in the code):

	- time to 1% creep strain
	- time to 2% creep strain
	- time to creep rupture