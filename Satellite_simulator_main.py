import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider
import warnings
import copy

# --- Physical Constants ---
MU_EARTH = 3.986004418e14
EARTH_RADIUS = 6378137.0
OMEGA_EARTH = 7.2921150e-5
C_LIGHT = 299792458.0

# --- Default Simulation Parameters ---
default_radius_init_m = EARTH_RADIUS + 600 * 1000
default_v_inertial_init = np.sqrt(MU_EARTH / default_radius_init_m)
default_velocity_init_ecef_y = default_v_inertial_init - OMEGA_EARTH * default_radius_init_m
default_initial_state_ecef_km = {
    'x_km': default_radius_init_m / 1000, 'y_km': 0.0, 'z_km': 0.0,
    'vx_mps': 0.0, 'vy_mps': default_velocity_init_ecef_y, 'vz_mps': 0.0
}
t_start_default = 0.0; t_end_default_min = 98.0; dt_default = 10.0
ue_beamwidth_default_deg = 30.0
visual_time_default_pct = 100.0

# --- Default UE Ground Point (Example: Mid-Latitude Europe ~45N, 10E) ---
default_lat_gp_deg = 45.0 # Changed Latitude
default_lon_gp_deg = 10.0 # Changed Longitude
default_lat_gp_rad = np.radians(default_lat_gp_deg)
default_lon_gp_rad = np.radians(default_lon_gp_deg)
default_r_gp_x = EARTH_RADIUS * np.cos(default_lat_gp_rad) * np.cos(default_lon_gp_rad)
default_r_gp_y = EARTH_RADIUS * np.cos(default_lat_gp_rad) * np.sin(default_lon_gp_rad)
default_r_gp_z = EARTH_RADIUS * np.sin(default_lat_gp_rad)
default_ground_point_ecef = [default_r_gp_x, default_r_gp_y, default_r_gp_z]
default_ground_point_ecef = [default_r_gp_x, 0, 0]

# --- User Input Section ---
print("--- Configuration ---")
use_defaults_sat = input("Use default initial satellite state? (y/n) [y]: ").lower()
initial_state_params = copy.deepcopy(default_initial_state_ecef_km)
if use_defaults_sat == 'n':
    print("\n--- Initial Satellite State (ECEF) ---")
    try:
        initial_state_params['x_km'] = float(input(f"  X (km) [def: {initial_state_params['x_km']:.1f}]: ") or initial_state_params['x_km'])
        initial_state_params['y_km'] = float(input(f"  Y (km) [def: {initial_state_params['y_km']:.1f}]: ") or initial_state_params['y_km'])
        initial_state_params['z_km'] = float(input(f"  Z (km) [def: {initial_state_params['z_km']:.1f}]: ") or initial_state_params['z_km'])
        initial_state_params['vx_mps'] = float(input(f"  Vx (m/s) [def: {initial_state_params['vx_mps']:.1f}]: ") or initial_state_params['vx_mps'])
        initial_state_params['vy_mps'] = float(input(f"  Vy (m/s) [def: {initial_state_params['vy_mps']:.1f}]: ") or initial_state_params['vy_mps'])
        initial_state_params['vz_mps'] = float(input(f"  Vz (m/s) [def: {initial_state_params['vz_mps']:.1f}]: ") or initial_state_params['vz_mps'])
    except ValueError: print("Invalid input. Using defaults for initial state."); initial_state_params = default_initial_state_ecef_km
else: print("Using default initial satellite state.")

while True:
    try:
        freq_ghz_str = input(f"Enter Carrier Frequency (GHz) [e.g., 20]: "); CARRIER_FREQUENCY_HZ = float(freq_ghz_str) * 1e9
        if CARRIER_FREQUENCY_HZ <= 0: raise ValueError("Frequency must be positive.")
        print(f"Using Carrier Frequency: {CARRIER_FREQUENCY_HZ/1e9:.2f} GHz"); break
    except ValueError as e: print(f"Invalid input: {e}. Please enter a positive number.")

ground_points_ecef_list = []
print("\n--- Define Ground Point(s) for Doppler/Beam Viz ---")
print(f"(Default GP #1 approx: Lat {default_lat_gp_deg:.1f}N, Lon {default_lon_gp_deg:.1f}E)")
print("(The FIRST point defined will be the origin for the UE beam visualization)")
use_default_gp = input(f"Use default ground point? (y/n) [y]: ").lower()
if use_default_gp != 'n':
    ground_points_ecef_list.append(default_ground_point_ecef)
    print(f"  Using default GP #1: [{default_ground_point_ecef[0]:.1f}, {default_ground_point_ecef[1]:.1f}, {default_ground_point_ecef[2]:.1f}]")
while True:
    prompt = "Add another fixed ground point? (y/n): " if ground_points_ecef_list else "Add the first fixed ground point? (y/n): "
    add_point = input(prompt).lower()
    if add_point != 'y':
        if not ground_points_ecef_list: print("Warning: Cannot visualize UE beam or calculate results without a ground point.")
        break
    print(f"\nEnter ECEF for ground point #{len(ground_points_ecef_list) + 1} (meters):")
    try:
        x_gp=float(input("  X (m): ")); y_gp=float(input("  Y (m): ")); z_gp=float(input("  Z (m): "))
        if np.linalg.norm([x_gp,y_gp,z_gp]) < EARTH_RADIUS*0.9: print("Warning: Point inside Earth.")
        ground_points_ecef_list.append([x_gp, y_gp, z_gp])
        print(f"  GP #{len(ground_points_ecef_list)} added: [{x_gp:.1f},{y_gp:.1f},{z_gp:.1f}]")
    except ValueError: print("Invalid input. Point not added.")
ground_points_ecef = np.array(ground_points_ecef_list); print("-" * 20)
omega_vec = np.array([0,0,OMEGA_EARTH]); v_grounds_ecef = [np.cross(omega_vec, r_g) for r_g in ground_points_ecef]
t_start = t_start_default; dt = dt_default

# --- Helper Functions ---
def ecef_to_eci(r_ecef, v_ecef, t):
    theta = OMEGA_EARTH * t; C, S = np.cos(theta), np.sin(theta)
    R_T = np.array([[C,-S,0],[S,C,0],[0,0,1]])
    r_eci = R_T @ r_ecef
    wXr = np.array([-OMEGA_EARTH*r_eci[1], OMEGA_EARTH*r_eci[0], 0])
    v_eci = R_T @ v_ecef + wXr
    return r_eci, v_eci

def eci_to_ecef(r_eci, v_eci, t):
    theta = OMEGA_EARTH * t; C, S = np.cos(theta), np.sin(theta)
    R = np.array([[C,S,0],[-S,C,0],[0,0,1]])
    r_ecef = R @ r_eci
    wXr = np.array([-OMEGA_EARTH*r_eci[1], OMEGA_EARTH*r_eci[0], 0])
    v_ecef = R @ (v_eci - wXr)
    return r_ecef, v_ecef

def derivatives(state_eci):
    r = state_eci[0:3]; r_mag = np.linalg.norm(r)
    if r_mag < EARTH_RADIUS * 0.9: return np.zeros(3)
    if r_mag < 1e-6: return np.zeros(3)
    return -MU_EARTH * r / (r_mag**3)

# --- Simulation Core ---
def calculate_orbit_with_vel(init_state_m,t_s,t_e_s,d_t):
    if t_e_s<=t_s:return np.array([]),np.array([]),np.array([])
    ts=np.arange(t_s,t_e_s+d_t,d_t)
    if len(ts)==0:return np.array([]),np.array([]),np.array([])
    pos,vel=[],[];r_ecf,v_ecf=init_state_m[0:3],init_state_m[3:6]
    for t_curr in ts:
        pos.append(r_ecf);vel.append(v_ecf)
        if np.linalg.norm(r_ecf)<EARTH_RADIUS*0.95:break
        r_eci,v_eci=ecef_to_eci(r_ecf,v_ecf,t_curr);acc=derivatives(np.concatenate((r_eci,v_eci)))
        r_eci_n=r_eci+v_eci*d_t;v_eci_n=v_eci+acc*d_t
        r_ecf,v_ecf=eci_to_ecef(r_eci_n,v_eci_n,t_curr+d_t)
    N=len(pos); return np.array(pos),np.array(vel),ts[:N]

# --- Geometry & Doppler ---
# calculate_projection_boundary removed as it's not used for UE beam or LOS plot

def calculate_doppler_and_los(r_sat, v_sat, r_ground, v_ground, fc):
    los_vector = r_ground - r_sat
    los_distance_m = np.linalg.norm(los_vector)
    relative_velocity_vector = v_sat - v_ground
    relative_speed_mps = np.linalg.norm(relative_velocity_vector)
    if los_distance_m < 1e-3:
        return 0.0, relative_speed_mps, 0.0
    los_unit_vector = los_vector / los_distance_m
    v_relative_los = np.dot(relative_velocity_vector, los_unit_vector)
    doppler_shift_hz = -(v_relative_los * fc) / C_LIGHT
    return doppler_shift_hz, relative_speed_mps, los_distance_m

def calculate_elevation_angle(r_sat, r_ground):
    vector_ground_to_sat = r_sat - r_ground
    vector_ground_to_zenith = r_ground
    dot_product = np.dot(vector_ground_to_sat, vector_ground_to_zenith)
    magnitude_product = np.linalg.norm(vector_ground_to_sat) * np.linalg.norm(vector_ground_to_zenith)
    if magnitude_product < 1e-9: return -90.0
    cos_zenith_angle = np.clip(dot_product / magnitude_product, -1.0, 1.0)
    zenith_angle_rad = np.arccos(cos_zenith_angle)
    elevation_rad = np.pi / 2.0 - zenith_angle_rad
    return np.degrees(elevation_rad)

def is_satellite_in_ue_beam(r_sat, r_ground, ue_beamwidth_deg):
    """Checks if satellite is within the UE's zenith-pointing beam cone."""
    vec_gp_zenith = r_ground
    mag_gp_zenith = np.linalg.norm(vec_gp_zenith)
    if mag_gp_zenith < 1e-9: return False
    unit_gp_zenith = vec_gp_zenith / mag_gp_zenith
    vec_gp_to_sat = r_sat - r_ground
    mag_gp_to_sat = np.linalg.norm(vec_gp_to_sat)
    if mag_gp_to_sat < 1e-9: return True
    unit_gp_to_sat = vec_gp_to_sat / mag_gp_to_sat
    dot_product = np.dot(unit_gp_zenith, unit_gp_to_sat)
    angle_rad = np.arccos(np.clip(dot_product, -1.0, 1.0))
    half_beamwidth_rad = np.radians(ue_beamwidth_deg) / 2.0
    return angle_rad <= half_beamwidth_rad

# --- Plot Setup ---
fig_3d = plt.figure("3D Orbit View", figsize=(9, 7)); ax_3d = fig_3d.add_axes([0.05,0.25,0.9,0.7],projection='3d')
u,v=np.mgrid[0:2*np.pi:100j,0:np.pi:50j]
ax_3d.plot_surface((EARTH_RADIUS/1000)*np.cos(u)*np.sin(v),(EARTH_RADIUS/1000)*np.sin(u)*np.sin(v),(EARTH_RADIUS/1000)*np.cos(v),color='b',alpha=0.3,rstride=1,cstride=1,zorder=1)
if ground_points_ecef.size>0: ax_3d.scatter(ground_points_ecef[:,0]/1000,ground_points_ecef[:,1]/1000,ground_points_ecef[:,2]/1000,c='m',marker='s',s=50,label='GPs',zorder=15)
orbit_line, = ax_3d.plot([],[],[],c='r',lw=1.5,label='Orbit',zorder=10)
# proj_line REMOVED
selected_sat_marker, = ax_3d.plot([],[],[],c='lime',marker='o',ms=8,label='Selected Sat Pos',zorder=11)
n_cone_lines=8; cone_lines=[ax_3d.plot([],[],[],c='c',ls='--',lw=0.8,zorder=4)[0] for _ in range(n_cone_lines)]

# --- Setup 2D Plots Figure (Now 6 subplots) ---
fig_2d, axs_2d = plt.subplots(6, 1, figsize=(9, 13), sharex=True) # Increased height for 6 plots
(ax_alt, ax_speed, ax_elev, ax_los, ax_delay, ax_dop) = axs_2d # Unpack axes, added ax_los, ax_delay
line_alt, = ax_alt.plot([], [], lw=1.5, label='Altitude')
line_speed, = ax_speed.plot([], [], lw=1.5, label='ECEF Speed')
lines_elev = [ax_elev.plot([], [], lw=1, marker='.', markersize=3, label=f'GP {i} Elev')[0] for i in range(len(ground_points_ecef))]
lines_los = [ax_los.plot([], [], lw=1, marker='.', markersize=3, label=f'GP {i} LOS')[0] for i in range(len(ground_points_ecef))] # Added LOS lines
lines_delay = [ax_delay.plot([], [], lw=1, marker='.', markersize=3, label=f'GP {i} Delay')[0] for i in range(len(ground_points_ecef))]
lines_dop = [ax_dop.plot([], [], lw=1, marker='.', markersize=3, label=f'GP {i} Doppler')[0] for i in range(len(ground_points_ecef))]

ax_alt.set_ylabel("Altitude (km)"); ax_alt.grid(True); ax_alt.legend(fontsize='small', loc='best')
ax_speed.set_ylabel("Speed (m/s)"); ax_speed.grid(True); ax_speed.legend(fontsize='small', loc='best')
ax_elev.set_ylabel("Elevation (°)"); ax_elev.grid(True); ax_elev.legend(fontsize='small', loc='best')
ax_los.set_ylabel("LOS Dist (km)"); ax_los.grid(True); ax_los.legend(fontsize='small', loc='best') # Added LOS axis
ax_delay.set_ylabel("Delay (ms)"); ax_delay.grid(True); ax_delay.legend(fontsize='small', loc='best')
ax_dop.set_ylabel("Doppler (kHz)"); ax_dop.set_xlabel("Time (min)"); ax_dop.grid(True); ax_dop.legend(fontsize='small', loc='best')
fig_2d.suptitle("Satellite Metrics vs. Time", fontsize=14)
fig_2d.tight_layout(rect=[0, 0.03, 1, 0.97]) # Adjust layout slightly more


# --- Global Caches ---
pos_cache, vel_cache, times_cache = np.array([]),np.array([]),np.array([])
last_initial_state, last_t_end_sec, last_ue_beam_deg = None,None,None
results_cache = {}

# --- Update Visuals (Zenith-Pointing UE Beam Cone) ---
def update_visuals_only(time_pct_val, ue_beam_deg_val):
    global cone_lines
    sat_vis, cone_vis = False, False

    for cl in cone_lines: cl.set_data_3d([], []); cl.set_visible(False)
    selected_sat_marker.set_data_3d([], [], []); selected_sat_marker.set_visible(False)
    # proj_line removed

    if pos_cache.size > 0 and ground_points_ecef.size > 0:
        gp_to_plot = ground_points_ecef[0]
        r_gp = gp_to_plot

        idx = min(max(0, int((time_pct_val / 100.0) * (len(pos_cache) - 1))), len(pos_cache) - 1)
        r_sat = pos_cache[idx]

        selected_sat_marker.set_data_3d([r_sat[0] / 1000], [r_sat[1] / 1000], [r_sat[2] / 1000])
        sat_vis = True
        selected_sat_marker.set_visible(sat_vis)

        if ue_beam_deg_val > 0:
            zenith_gp_unit = r_gp / np.linalg.norm(r_gp)
            z_axis = np.array([0,0,1.0]); U_gp = np.cross(z_axis, zenith_gp_unit)
            if np.linalg.norm(U_gp) < 1e-9: U_gp = np.cross([0,1,0], zenith_gp_unit)
            if np.linalg.norm(U_gp) < 1e-9: U_gp = np.cross([1,0,0], zenith_gp_unit)
            U_gp /= np.linalg.norm(U_gp)
            V_gp = np.cross(zenith_gp_unit, U_gp)

            beta_ue = np.radians(ue_beam_deg_val) / 2.0
            edge_rays_gp = []
            for i in range(n_cone_lines):
                alpha = 2 * np.pi * i / n_cone_lines
                term_uv = np.cos(alpha) * U_gp + np.sin(alpha) * V_gp
                V_ray = np.cos(beta_ue) * zenith_gp_unit + np.sin(beta_ue) * term_uv
                edge_rays_gp.append(V_ray)

            gp_pos_km = r_gp / 1000.0
            cone_draw_length_m = 1500 * 1000
            cone_vis = True
            for i in range(n_cone_lines):
                ray_endpoint_m = r_gp + cone_draw_length_m * edge_rays_gp[i]
                ray_endpoint_km = ray_endpoint_m / 1000.0
                cone_lines[i].set_data_3d([gp_pos_km[0], ray_endpoint_km[0]],
                                          [gp_pos_km[1], ray_endpoint_km[1]],
                                          [gp_pos_km[2], ray_endpoint_km[2]])
                cone_lines[i].set_visible(True)

    if not cone_vis:
        for cl in cone_lines: cl.set_visible(False)

# --- Print Results Table ---
def print_results_table(res_cache):
    if not res_cache:
        if ground_points_ecef.size>0: print("\n--- Results ---\nNo satellite within UE beam during simulation.\n"+"-"*100)
        return
    print("\n--- Results Table ---")
    hdr="Time(min)|SatSpd(m/s)|RelVel(m/s)|GP|Elev(°)|SatAlt(km)|LOS(km) |Delay(ms)|Doppler(kHz)"
    sep="-"*len(hdr)
    all_res=[item+(gp_idx,) for gp_idx,data in res_cache.items() for item in data]
    all_res.sort(key=lambda x:x[0])
    if all_res:
        print(hdr); print(sep)
        for t, dop, el, s_spd, s_alt, r_vel, los_km, delay_ms, gp in all_res:
             print(f"{t/60:9.2f}|{s_spd:12.1f}|{r_vel:11.1f}|{gp:2d}|{el:7.1f}|{s_alt:10.1f}|{los_km:8.0f}|{delay_ms:9.2f}|{dop/1000:12.3f}")
    elif ground_points_ecef.size>0: print("No satellite within UE beam during simulation.")
    print(sep)

# --- Main Update Function ---
def update(val=None):
    global pos_cache, vel_cache, times_cache
    global last_initial_state, last_t_end_sec, last_ue_beam_deg, results_cache

    curr_x_km = s_x.val; curr_y_km = s_y.val; curr_z_km = s_z.val
    curr_vx_mps = s_vx.val; curr_vy_mps = s_vy.val; curr_vz_mps = s_vz.val
    curr_t_end_min = s_tend.val
    curr_ue_beam_deg = s_ue_beam.val
    curr_time_pct = s_time_pct.val

    current_init_state_m = np.array([
        curr_x_km*1000, curr_y_km*1000, curr_z_km*1000,
        curr_vx_mps, curr_vy_mps, curr_vz_mps
    ])
    current_t_end_s = curr_t_end_min * 60

    orbit_params_changed = (last_initial_state is None or
                            not np.array_equal(current_init_state_m, last_initial_state) or
                            current_t_end_s != last_t_end_sec)
    ue_beam_param_changed = (last_ue_beam_deg is None or curr_ue_beam_deg != last_ue_beam_deg)
    recalculate_results = orbit_params_changed or ue_beam_param_changed

    if orbit_params_changed:
        print("\nRecalculating orbit...")
        pos_cache, vel_cache, times_cache = calculate_orbit_with_vel(current_init_state_m, t_start, current_t_end_s, dt)
        last_initial_state = copy.deepcopy(current_init_state_m)
        last_t_end_sec = current_t_end_s
        print("Orbit calculation complete.")

        orbit_vis = pos_cache.size > 0
        if orbit_vis: orbit_line.set_data_3d(pos_cache[:,0]/1000, pos_cache[:,1]/1000, pos_cache[:,2]/1000)
        else: orbit_line.set_data_3d([],[],[])
        orbit_line.set_visible(orbit_vis)

        if orbit_vis:
            try:
                max_r = np.max(np.linalg.norm(pos_cache/1000,axis=1))*1.1
                lim = max(max_r, plot_radius_default_init * 0.5)
                ax_3d.set_xlim([-lim,lim]);ax_3d.set_ylim([-lim,lim]);ax_3d.set_zlim([-lim,lim])
            except ValueError: pass

    if recalculate_results:
        print("Calculating Doppler/Elevation/etc (based on UE Beamwidth)...")
        results_cache = {i:[] for i in range(len(ground_points_ecef))}
        if ground_points_ecef.size>0 and pos_cache.size>0:
            for i,t in enumerate(times_cache):
                r_s,v_s = pos_cache[i],vel_cache[i]
                s_spd=np.linalg.norm(v_s); s_alt=(np.linalg.norm(r_s)-EARTH_RADIUS)/1000
                for gp_i,r_g in enumerate(ground_points_ecef):
                    v_g=v_grounds_ecef[gp_i]
                    if is_satellite_in_ue_beam(r_s,r_g,curr_ue_beam_deg):
                        dop_hz,r_vel,los_m = calculate_doppler_and_los(r_s,v_s,r_g,v_g,CARRIER_FREQUENCY_HZ)
                        el_deg = calculate_elevation_angle(r_s,r_g)
                        delay_ms = (los_m / C_LIGHT) * 1000.0
                        results_cache[gp_i].append((t,dop_hz,el_deg,s_spd,s_alt,r_vel, los_m / 1000.0, delay_ms))
        print("Calculations complete.")
        print_results_table(results_cache)
        last_ue_beam_deg = curr_ue_beam_deg

        # --- Update 2D Plots ---
        times_min_cache = times_cache / 60.0 if times_cache.size > 0 else np.array([])
        if pos_cache.size > 0: line_alt.set_data(times_min_cache, (np.linalg.norm(pos_cache, axis=1) - EARTH_RADIUS) / 1000.0); ax_alt.relim(); ax_alt.autoscale_view()
        else: line_alt.set_data([],[])
        if vel_cache.size > 0: line_speed.set_data(times_min_cache, np.linalg.norm(vel_cache, axis=1)); ax_speed.relim(); ax_speed.autoscale_view()
        else: line_speed.set_data([],[])

        min_elev_overall,max_elev_overall=-90,90; min_los_overall,max_los_overall=0,0
        min_delay_overall,max_delay_overall=0,0; min_dop_overall,max_dop_overall = 0,0
        data_plotted_elev, data_plotted_dop, data_plotted_delay, data_plotted_los = False, False, False, False

        for gp_i in range(len(ground_points_ecef)):
            if gp_i in results_cache and results_cache[gp_i]:
                # Unpack the full tuple including LOS and delay
                gp_times, gp_dops, gp_els, _, _, _, gp_los_km, gp_delays_ms = zip(*results_cache[gp_i])
                gp_times_min = np.array(gp_times)/60.0; gp_dops_khz = np.array(gp_dops)/1000.0

                lines_elev[gp_i].set_data(gp_times_min, gp_els)
                lines_dop[gp_i].set_data(gp_times_min, gp_dops_khz)
                lines_delay[gp_i].set_data(gp_times_min, gp_delays_ms)
                lines_los[gp_i].set_data(gp_times_min, gp_los_km) # Set LOS data

                min_e, max_e = np.min(gp_els), np.max(gp_els); min_elev_overall=min(min_elev_overall,min_e) if data_plotted_elev else min_e; max_elev_overall=max(max_elev_overall,max_e) if data_plotted_elev else max_e; data_plotted_elev=True
                min_d, max_d = np.min(gp_dops_khz), np.max(gp_dops_khz); min_dop_overall=min(min_dop_overall,min_d) if data_plotted_dop else min_d; max_dop_overall=max(max_dop_overall,max_d) if data_plotted_dop else max_d; data_plotted_dop=True
                min_del, max_del = np.min(gp_delays_ms), np.max(gp_delays_ms); min_delay_overall=min(min_delay_overall,min_del) if data_plotted_delay else min_del; max_delay_overall=max(max_delay_overall,max_del) if data_plotted_delay else max_del; data_plotted_delay=True
                min_l, max_l = np.min(gp_los_km), np.max(gp_los_km); min_los_overall=min(min_los_overall,min_l) if data_plotted_los else min_l; max_los_overall=max(max_los_overall,max_l) if data_plotted_los else max_l; data_plotted_los=True

            else:
                lines_elev[gp_i].set_data([],[]); lines_dop[gp_i].set_data([],[]); lines_delay[gp_i].set_data([],[]); lines_los[gp_i].set_data([],[]) # Clear LOS line too

        if data_plotted_elev: margin_e=(max_elev_overall-min_elev_overall)*0.1 if (max_elev_overall-min_elev_overall)>1e-3 else 5; ax_elev.set_ylim([min(-5,min_elev_overall-margin_e),max(95,max_elev_overall+margin_e)]); ax_elev.legend(fontsize='small', loc='best')
        else: ax_elev.set_ylim([-10,95])
        if data_plotted_dop: margin_d=(max_dop_overall-min_dop_overall)*0.1 if (max_dop_overall-min_dop_overall)>1e-3 else 0.5; ax_dop.set_ylim([min_dop_overall-margin_d-0.5,max_dop_overall+margin_d+0.5]); ax_dop.legend(fontsize='small', loc='best')
        else: ax_dop.set_ylim([-1,1])
        if data_plotted_delay: margin_del=(max_delay_overall-min_delay_overall)*0.1 if (max_delay_overall-min_delay_overall)>0.1 else 1; ax_delay.set_ylim([max(0,min_delay_overall-margin_del),max_delay_overall+margin_del]); ax_delay.legend(fontsize='small', loc='best')
        else: ax_delay.set_ylim([0,10])
        if data_plotted_los: margin_l=(max_los_overall-min_los_overall)*0.05; ax_los.set_ylim([max(0,min_los_overall-margin_l),max_los_overall+margin_l]); ax_los.legend(fontsize='small', loc='best')
        else: ax_los.set_ylim([0, 5000]) # Default LOS range
        fig_2d.canvas.draw_idle()

    if pos_cache.size > 0: update_visuals_only(curr_time_pct, curr_ue_beam_deg)
    else: update_visuals_only(0, 0)

    ax_3d.set_title(f'Orbit & UE Zenith Beam (UE BW:{curr_ue_beam_deg:.1f}°, Time:{curr_time_pct:.1f}%)')
    fig_3d.canvas.draw_idle()

# Plot limits & labels
plot_radius_default_init = np.linalg.norm([initial_state_params['x_km'], initial_state_params['y_km'], initial_state_params['z_km']]) * 1.5
ax_3d.set_xlim([-plot_radius_default_init, plot_radius_default_init]); ax_3d.set_ylim([-plot_radius_default_init, plot_radius_default_init]); ax_3d.set_zlim([-plot_radius_default_init, plot_radius_default_init])
ax_3d.set_xlabel('X(km)'); ax_3d.set_ylabel('Y(km)'); ax_3d.set_zlabel('Z(km)')
ax_3d.set_title('Satellite Orbit'); ax_3d.set_aspect('equal',adjustable='box'); ax_3d.legend(loc='upper right',fontsize='small')

# Sliders
axcolor='lightgoldenrodyellow'; slider_h=0.020; slider_vspace=0.006
num_sliders = 9 # X,Y,Z, Vx,Vy,Vz, Tend, UEBeam, TimePct
slider_start_y = 0.01
slider_axes=[fig_3d.add_axes([0.15,slider_start_y+(num_sliders-1-i)*(slider_h+slider_vspace),0.7,slider_h],fc=axcolor) for i in range(num_sliders)]
s_x_ax, s_y_ax, s_z_ax, s_vx_ax, s_vy_ax, s_vz_ax, s_tend_ax, s_ue_beam_ax, s_time_pct_ax = slider_axes

# Define Slider Ranges
p_range=plot_radius_default_init; p_min,p_max=-p_range,p_range
v_range=10000; v_min,v_max=-v_range,v_range
t_min,t_max=5.0,200.0; tp_min,tp_max=0.0,100.0
ue_b_min, ue_b_max = 0.1, 160.0

# Create Sliders
s_x = Slider(s_x_ax,'X(km)',p_min,p_max,valinit=initial_state_params['x_km'],valstep=10)
s_y = Slider(s_y_ax,'Y(km)',p_min,p_max,valinit=initial_state_params['y_km'],valstep=10)
s_z = Slider(s_z_ax,'Z(km)',p_min,p_max,valinit=initial_state_params['z_km'],valstep=10)
s_vx = Slider(s_vx_ax,'Vx(m/s)',v_min,v_max,valinit=initial_state_params['vx_mps'],valstep=10)
s_vy = Slider(s_vy_ax,'Vy(m/s)',v_min,v_max,valinit=initial_state_params['vy_mps'],valstep=10)
s_vz = Slider(s_vz_ax,'Vz(m/s)',v_min,v_max,valinit=initial_state_params['vz_mps'],valstep=10)
s_tend = Slider(s_tend_ax,'End Time(min)',t_min,t_max,valinit=t_end_default_min,valstep=1)
s_ue_beam = Slider(s_ue_beam_ax,'UE BW (° - Doppler/Visual)',ue_b_min, ue_b_max,valinit=ue_beamwidth_default_deg,valstep=1.0)
s_time_pct = Slider(s_time_pct_ax,'Visual Time(%)',tp_min,tp_max,valinit=visual_time_default_pct,valstep=0.5)

# Initial run
update()

# Connect sliders
s_x.on_changed(update); s_y.on_changed(update); s_z.on_changed(update)
s_vx.on_changed(update); s_vy.on_changed(update); s_vz.on_changed(update)
s_tend.on_changed(update); s_ue_beam.on_changed(update); s_time_pct.on_changed(update)

plt.show()