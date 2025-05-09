________________________________________
Interactive Satellite Orbit Simulation with Beam Projection and Doppler Calculation
Version: Based on the Python script developed on 2024-07-27
________________________________________
1. Introduction
This document describes a Python script that simulates the orbit of a Low Earth Orbit (LEO) satellite around the Earth. It provides an interactive 3D visualization using Matplotlib, allowing users to adjust the satellite's initial state (position and velocity), simulation duration, and antenna beamwidth via sliders.
The simulation calculates and displays:
•	The satellite's trajectory in the Earth-Cantered, Earth-Fixed (ECEF) coordinate frame.
•	The projected boundary of the satellite's antenna beam cone onto the Earth's surface at a user-selected point in time.
•	A visual representation of the beam cone itself originating from the satellite.
•	The Doppler shift experienced by fixed ground points on Earth when they are within the satellite's beam cone.
The primary purpose of this tool is educational, demonstrating the interplay between orbital mechanics, coordinate systems, beam geometry, and the Doppler effect in the context of satellite communications, particularly relevant to Non-Terrestrial Networks (NTN).
________________________________________
2. Core Concepts
Several fundamental concepts underpin this simulation:
•	Orbital Mechanics (Two-Body Problem): The simulation simplifies the satellite's motion by considering only the gravitational force between the Earth and the satellite (the two-body problem). It ignores perturbations like atmospheric drag, gravitational anomalies from Earth's non-uniformity (e.g., J2 effect), solar radiation pressure, and third-body influences (Moon, Sun). The resulting motion under ideal conditions follows Kepler's laws.
•	Coordinate Systems:
o	Earth-Centered Inertial (ECI): A non-rotating frame with its origin at the Earth's center. Its axes are typically fixed relative to distant stars (e.g., pointing towards the vernal equinox). Physics laws (like gravity and motion) are simplest to apply in this frame.
o	Earth-Centered, Earth-Fixed (ECEF): A frame with its origin at the Earth's center that rotates along with the Earth. Its X-axis typically points towards the intersection of the Prime Meridian and the Equator, the Z-axis points towards the North Pole, and the Y-axis completes the right-handed system. Locations on Earth (like ground stations) are fixed in this frame.
o	Transformation: Moving between these frames requires accounting for the Earth's rotation (OMEGA_EARTH).
•	Beam Geometry:
o	Nadir Pointing: The simulation assumes the satellite's antenna beam is always pointed directly towards the center of the Earth (the nadir direction).
o	Cone Model: The beam is modeled as a simple cone defined by its beamwidth (typically the angle where power drops by 3dB).
o	Projection: The intersection of this cone with the spherical Earth surface defines the projected beam footprint boundary, which is generally an ellipse.
•	Doppler Effect: When there is relative motion between a transmitter (satellite) and a receiver (ground point) along the line connecting them (line-of-sight), the received frequency is shifted relative to the transmitted frequency. The magnitude and sign of the shift depend on the relative velocity along the line-of-sight.
________________________________________
3. Code Structure Breakdown
The Python script is organized as follows:
•	Import Statements: Imports necessary libraries (numpy, matplotlib, Slider, warnings, copy).
•	Physical Constants: Defines key physical values (MU_EARTH, EARTH_RADIUS, OMEGA_EARTH, C_LIGHT).
•	Simulation Parameters: Sets default values for initial state, time steps, and beamwidth.
•	User Input: Prompts the user to enter the carrier frequency and the ECEF coordinates of fixed ground points. Pre-calculates the velocity of these ground points due to Earth's rotation.
•	Helper Functions:
o	ecef_to_eci(r_ecef, v_ecef, t): Transforms position and velocity from ECEF to ECI at a given time t.
o	eci_to_ecef(r_eci, v_eci, t): Transforms position and velocity from ECI to ECEF at a given time t.
o	derivatives(state_eci): Calculates the gravitational acceleration vector acting on the satellite in the ECI frame.
•	Simulation Core:
o	calculate_orbit_with_vel(...): Propagates the satellite's state vector over time using numerical integration. It takes the initial ECEF state, start/end times, and time step (dt). It returns arrays containing the ECEF positions, ECEF velocities, and corresponding times for each step.
•	Geometry & Doppler Calculations:
o	calculate_projection_boundary(...): Calculates the ECEF coordinates of the points forming the boundary where the beam cone intersects the Earth sphere. It also returns the direction vectors of the rays forming the cone edge and the maximum distance from the satellite to a projection point (used for drawing the cone).
o	is_ground_point_in_beam(...): Checks if a given ground point lies geometrically within the beam cone originating from the satellite.
o	calculate_doppler_at_point(...): Calculates the Doppler shift in Hz for a specific satellite and ground point state, given the carrier frequency.
•	Plotting Setup:
o	Creates the Matplotlib figure and 3D axes.
o	Plots the Earth as a blue sphere.
o	Plots the fixed ground points entered by the user (magenta squares).
o	Initializes empty plot objects (orbit_line, projection_line, cone_lines) that will be updated dynamically.
•	Global Variables: Used to cache the results of the last orbit calculation (ecef_positions_cache, etc.) to avoid redundant computations when only display parameters change.
•	Update Logic:
o	draw_plots(selected_sat_pos, beam_deg): Updates the visual elements (projection boundary line, cone lines) based on a specific satellite position and beamwidth. This function does not recalculate the orbit.
o	update(val): The main callback function triggered by any slider change.
	Reads all slider values.
	Determines if the orbit needs recalculation based on changes to initial state or end time sliders.
	If recalculation is needed:
	Calls calculate_orbit_with_vel.
	Updates the cache.
	Updates the orbit_line plot.
	Performs the Doppler calculation loop for the entire new trajectory and prints results to the console.
	Selects the satellite position from the cache based on the "Visual Time (%)" slider.
	Calls draw_plots to update the projection and cone visuals for the selected time.
	Updates the plot title.
	Redraws the figure canvas.
•	Slider Setup: Creates the slider axes and the Slider widgets, linking them to the update function.
•	Initial Draw: Calls update(None) once to perform the initial orbit calculation, Doppler print, and plot setup before displaying the figure.
•	plt.show(): Displays the interactive plot window.
________________________________________
4. Physics and Equations
•	Gravitational Acceleration (ECI Frame):
Based on Newton's Law of Universal Gravitation, the acceleration a of the satellite due to Earth's gravity is:
a = - (GM / |r|^2) * (r / |r|) = - μ * r / |r|^3
Where:
* μ = GM is Earth's standard gravitational parameter (MU_EARTH).
* r is the position vector of the satellite relative to Earth's center in the ECI frame.
* |r| is the magnitude of the position vector (distance from Earth's center).
Reference: Any standard physics or astrodynamics textbook (e.g., "Fundamentals of Astrodynamics" by Bate, Mueller, White; "Orbital Mechanics for Engineering Students" by Howard D. Curtis).
(Implemented in derivatives function)
•	Numerical Integration (Euler Method):
To find the position and velocity at the next time step (t + dt), the simple Euler method is used:
v_new = v_old + a_old * dt
r_new = r_old + v_old * dt
Where:
* v_old, r_old are velocity and position at time t.
* a_old is the acceleration calculated using r_old (and potentially v_old if non-gravitational forces were included).
* dt is the time step.
This is a first-order method; higher-order methods (like Runge-Kutta 4) provide better accuracy for the same dt but are more complex to implement.
Reference: Numerical analysis textbooks.
(Implemented in calculate_orbit_with_vel function)
•	Coordinate Transformation (ECI <=> ECEF):
The transformation involves a rotation around the Earth's Z-axis (polar axis) by an angle θ = ωe * t.
o	Rotation Matrix (ECI to ECEF):
o	      R(θ) = |  cos(θ)   sin(θ)   0 |
o	       | -sin(θ)   cos(θ)   0 |
o	       |    0        0      1 |
    
o	Position Transformation:
r_ecef = R(θ) * r_eci
r_eci = R(θ)^T * r_ecef (where R^T is the transpose)
o	Velocity Transformation: This requires accounting for the rotating frame. The relationship is:
v_eci = R(θ)^T * v_ecef + ω × r_eci
v_ecef = R(θ) * (v_eci - ω × r_eci)
Where ω = [0, 0, ωe] is the Earth's angular velocity vector in the ECI frame. The ω × r_eci term accounts for the velocity induced by the frame's rotation.
Reference: Astrodynamics textbooks (e.g., Vallado, Curtis), texts on classical mechanics.
(Implemented in ecef_to_eci and eci_to_ecef functions)
•	Ray-Sphere Intersection:
To find where a ray starting from the satellite (S) in direction (V) intersects the Earth sphere (|P|^2 = R_e^2), we solve for the distance d along the ray:
P = S + d * V
|S + d * V|^2 = R_e^2
Expanding this leads to a quadratic equation in d:
|V|^2 * d^2 + 2 * (S ⋅ V) * d + (|S|^2 - R_e^2) = 0
Since V is a unit vector (|V|^2 = 1), this simplifies to:
d^2 + (2 * S ⋅ V) * d + (|S|^2 - R_e^2) = 0
This is Ad^2 + Bd + C = 0 where A=1, B = 2 * S ⋅ V, C = |S|^2 - R_e^2. The code solves this using the quadratic formula d = (-B ± sqrt(B^2 - 4AC)) / 2A. The smallest positive real root d corresponds to the intersection point closest to the satellite.
Reference: Computer graphics textbooks, analytic geometry.
(Implemented in calculate_projection_boundary function)
•	Beam Cone Angle Check:
To determine if a ground point G is within the beam cone originating from satellite S with half-beamwidth β and pointing towards nadir (N_sat), we check the angle α between the satellite-to-nadir vector (-S) and the satellite-to-ground vector (G - S):
α = arccos( ((-S) ⋅ (G - S)) / (|-S| * |G - S|) )
The point is inside the beam if α <= β.
Reference: Basic vector geometry.
(Implemented in is_ground_point_in_beam function)
•	Doppler Shift:
The first-order Doppler shift fd is calculated as:
fd = - v_rel_los / λ = - (v_rel ⋅ u_los) * fc / c
Where:
* v_rel = v_sat - v_ground is the relative velocity vector between satellite and ground point (in the same frame, ECEF).
* u_los = (r_ground - r_sat) / |r_ground - r_sat| is the unit vector along the line-of-sight from satellite to ground point.
* v_rel_los = v_rel ⋅ u_los is the component of relative velocity along the line-of-sight (positive means increasing distance).
* λ = c / fc is the wavelength of the carrier signal.
* fc is the carrier frequency (CARRIER_FREQUENCY_HZ).
* c is the speed of light (C_LIGHT).
The negative sign indicates that frequency increases (positive Doppler) when the distance is decreasing (negative v_rel_los).
Reference: Physics textbooks (wave phenomena), radio communication engineering textbooks.
(Implemented in calculate_doppler_at_point function)
________________________________________
5. How to Run and Use
1.	Save: Save the code as a Python file (e.g., satellite_sim.py).
2.	Install Dependencies: Ensure you have numpy and matplotlib installed (pip install numpy matplotlib).
3.	Run: Execute the script from your terminal (python satellite_sim.py).
4.	Input: The script will first prompt you to enter:
o	The carrier frequency in GHz.
o	Whether to add fixed ground points (y/n).
o	If 'y', the ECEF X, Y, Z coordinates (in meters) for each point.
5.	Interact:
o	A 3D plot window will appear showing the Earth, ground points (if any), and the initial calculated satellite orbit.
o	Use the sliders at the bottom to adjust:
	Initial ECEF position (X, Y, Z in km) and velocity (Vx, Vy, Vz in m/s).
	Simulation duration ("End Time" in minutes).
	Antenna beamwidth in degrees.
	The time point ("Visual Time (%)") for displaying the beam projection and cone.
o	When sliders controlling the initial state or end time are changed, the orbit is recalculated, the orbit plot updates, and the Doppler shift calculations are run again and printed to the console for the entire new trajectory.
o	When sliders controlling beamwidth or visual time are changed, only the visual elements (orange projection line, cyan cone lines) are updated based on the cached orbit data. The Doppler calculation is not rerun.
o	Use the standard Matplotlib tools to rotate, pan, and zoom the 3D plot for better visualization.
6.	Output: Doppler shift values (in kHz) are printed to the console whenever a ground point is within the beam cone during the simulation pass associated with the last orbit recalculation.
________________________________________
6. Limitations and Potential Improvements
•	Simplified Physics: Uses a basic two-body gravity model. Adding perturbations (J2, drag, etc.) would increase realism but also complexity.
•	Euler Integrator: Prone to accumulating errors over long simulations. Using a higher-order integrator (e.g., RK4 via scipy.integrate.solve_ivp) is recommended for better accuracy.
•	Nadir Pointing: Assumes the beam is always pointing directly down. Real systems have specific pointing strategies and potential pointing errors.
•	Spherical Earth: Uses a spherical Earth model for some geometry calculations (like the Earth plot and ray intersection). A more accurate ellipsoidal model (WGS84) could be used.
•	Instantaneous Doppler: Calculates Doppler at discrete time steps. It doesn't model the continuous change or effects like Doppler rate.
•	Atmospheric Effects: Ignores signal delay, refraction, or attenuation caused by the atmosphere.
•	Visualization: The cone is represented by lines; a transparent surface could be used but might impact performance. The projection is only shown at one time instant.
•	Performance: Recalculating the orbit can be slow for long durations or small time steps.
________________________________________
7. Conclusion
This interactive simulation provides a valuable tool for visualizing LEO satellite orbits in the ECEF frame and understanding fundamental concepts like beam projection and the Doppler effect. While based on simplified models, it allows users to explore how changes in initial conditions, orbit duration, and beamwidth affect the satellite's path, its ground coverage footprint, and the resulting frequency shift observed at specific ground locations. It serves as a foundation upon which more complex and realistic NTN simulations can be built.

