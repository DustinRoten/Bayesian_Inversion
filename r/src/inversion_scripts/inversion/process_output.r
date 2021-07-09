#This program is in charge of running the necessary scripts that 
#post processes output from our inversion code, and saves the data 
#either as a netcdf or as an rds file, depending on whether the 
#output is gridded or not. Time series and scatter plots are then
#saved in ./analysis_wrap/plots/ while rds and netcdf files are 
#saved to the ./out directory. 
#
#DVM 3/8/2020

#1) Write our priopr emissions to a netcdf file
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("Writing prior emissons to netcdf")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("./process_output/write_src/write_sprior.r")

#2) Write our posterior emissions to a netcdf file
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("Writing prior emissons to netcdf")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("./process_output/write_src/write_shat.r")

#3) Write our prior emission uncertainty to a netcdf file
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("Writing prior emission uncertainty to a netcdf file")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("./process_output/write_src/write_sigma.r")

#4. Write uncertainty reduction to a netcdf file
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("Write percent reduction in uncertainty")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("./process_output/write_src/write_percent_uncertainty_reduction.r")

#5. Write average STILT footprint to a netcdf file
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("Write average STILT footprint")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("./process_output/write_src/write_footprint.r")

#6. Plot timeseries of emissions
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("Plotting timeseries of emissions")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("./process_output/plot_src/plot_flux_timeseries.r")

#7. Convolve receptors
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("Convolving receptors")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("./process_output/write_src/convolve_receptor.r")

#8. Plot TRAX data and save data as kmz file
trax.line = "red_line"
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("Plot TRAX red data and save to kmz")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("./process_output/plot_src/trax_analysis.r")

trax.line = "green_line"
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("Plot TRAX green data and save to kmz")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("./process_output/plot_src/trax_analysis.r")

#9. Plot scatter plot
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("Plot Scatterplot")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("./process_output/plot_src/plot_scatt_plot.r")

#10. Plot diurnal cycle
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("Plot diurnal cycle")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("./process_output/plot_src/plot_diurnal_cycle.r")

