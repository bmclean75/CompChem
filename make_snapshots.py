
#Tested on pvpython version 5.11.1

import paraview.simple as pvs

import os

num_files = 500
min_e = 1.6
max_e = 3.2
min_model_dimensions = 500.0e-12

# Get the current working directory
current_working_directory = os.getcwd()

# Specify relative paths
data_relative_path = "output/"  # Adjust as needed
output_relative_path = "images/"  # Adjust as needed

# Construct the full paths by joining them with the current working directory
#data_path = current_working_directory + data_relative_path #os.path.join(current_working_directory, data_relative_path)
data_path = os.path.join(current_working_directory, data_relative_path)
output_path = os.path.join(current_working_directory, output_relative_path)

# Ensure the output directory exists or create it
if not os.path.exists(output_path):
	os.makedirs(output_path)

# Function to load a CSV file and create a visualization
def visualize_csv(filename, timestep):
	# Load the CSV file
	csv_reader = pvs.CSVReader(FileName=[filename])
	csv_reader.HaveHeaders = True
	
	# Create a table to points filter
	table_to_points = pvs.TableToPoints(Input=csv_reader)
	table_to_points.XColumn = "x"
	table_to_points.YColumn = "y"
	table_to_points.ZColumn = "z"
	
	# Apply a Glyph filter to represent atoms as spheres
	glyph_filter = pvs.Glyph(Input=table_to_points,
							GlyphType='Sphere')
	glyph_filter.GlyphType.Radius = 1  # Default radius, will be scaled by 'relsize'
	glyph_filter.GlyphMode = 'All Points'

	# Increase the resolution of the sphere glyphs for smoother appearance
	glyph_filter.GlyphType.ThetaResolution = 64  # Default is usually 8 or 16
	glyph_filter.GlyphType.PhiResolution = 64  # Default is usually 8 or 16

	
	# Scale spheres based on the 'relsize' attribute from the CSV
	glyph_filter.ScaleArray = ['POINTS', 'relsize']
	glyph_filter.ScaleFactor = min_model_dimensions / 150  # This factor will be multiplied by each point's 'relsize' value
	
	# Setup color mapping for 'e' values
	lookup_table = pvs.GetColorTransferFunction('e')
	lookup_table.RGBPoints = [min_e, 0.0, 0.0, 1.0, max_e, 1.0, 0.0, 0.0]  # Example, adjust according to your data's range
	lookup_table.ColorSpace = 'RGB'

	# Get active view or create one if not present
	render_view = pvs.GetActiveViewOrCreate('RenderView')

	# Set the render view size for higher resolution output
	render_view.ViewSize = [1080, 1080]  # For example, to set the resolution to 1920x1080 pixels

	# Display the glyphs in the render view
	glyph_display = pvs.Show(glyph_filter, render_view)

	# Directly apply the lookup table to the glyphs' representation
	glyph_display.LookupTable = lookup_table
	glyph_display.ColorArrayName = ['POINTS', 'e']

	# Automatically adjust the camera to show all data and get initial camera settings
	pvs.ResetCamera(render_view)
	initial_camera_position = render_view.CameraPosition
	initial_focal_point = render_view.CameraFocalPoint
	initial_view_up = render_view.CameraViewUp

	# Calculate new camera position for zooming
	zoom_factor = 1.5  # Example zoom factor; adjust as needed
	new_camera_position = [(initial_focal_point[i] + (initial_camera_position[i] - initial_focal_point[i]) / zoom_factor) for i in range(3)]

	# Apply the new camera position
	render_view.CameraPosition = new_camera_position
	render_view.CameraFocalPoint = initial_focal_point
	render_view.CameraViewUp = initial_view_up

	# Update the view to apply the new camera settings
	pvs.Render()

	# Save screenshot for the current timestep
	#screenshot_filename = f"{output_path}timestep_{timestep:04d}.png"
	screenshot_filename = f"{output_path}timestep_{timestep}.png"
	pvs.SaveScreenshot(screenshot_filename, render_view)
	print(f"Saved: {screenshot_filename}")

	# Delete the current glyph filter after rendering to keep the pipeline clean
	pvs.Delete(glyph_filter)
	pvs.Delete(table_to_points)
	pvs.Delete(csv_reader)

	# Ensure the deletion is reflected in the view
	pvs.Render()

# Main loop to process each file
for i in range(num_files):
	filename = f"{data_path}{i}.csv"
	visualize_csv(filename, i)

print("Visualization complete.")
