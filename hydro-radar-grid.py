from turtle import color
from matplotlib import markers
import matplotlib.pyplot as plt
import numpy as np

PIVOT_HEIGHT = 0.56 + 3.0 # m. Vertical distance between pivot and ground level.
ANTENNA_CENTER_POSITION = 0.12 # m. Distance between pivot and center of RX antenna.
MAX_SCAN_ANGLE = np.deg2rad(15) # rad
INITIAL_TILT_ANGLE = np.deg2rad(10) # rad
FINAL_TILT_ANGLE = np.deg2rad(40) # rad
TILT_POSITIONS = 4
SCANNING_DIRECTIONS = 5
BEAMWIDTH_AZIMUTH = np.deg2rad(6) # rad
BEAMWIDTH_ELEVATION = np.deg2rad(9) # rad

# Only for my particular radar prototype
#    /
#   /. RX antenna
#  /    .
# *pivot   .
# |           .
# |              . d (@ beam direction)
# |                 .
# |                    .
# |                       .
# ____                       ._____ground level
#    <----------------------->
#                l
antennaHeight = PIVOT_HEIGHT + ANTENNA_CENTER_POSITION * np.cos(INITIAL_TILT_ANGLE)
print(f"Radar height: {antennaHeight:.1f}")
tiltAngles = np.linspace(start=INITIAL_TILT_ANGLE, stop=FINAL_TILT_ANGLE, num=TILT_POSITIONS, endpoint=True)
directions = np.linspace(start=MAX_SCAN_ANGLE, stop=-MAX_SCAN_ANGLE, num=SCANNING_DIRECTIONS, endpoint=True)
grid_parameters = np.zeros((TILT_POSITIONS, SCANNING_DIRECTIONS, 8))
# For each pair (tiltAngle, direction) the following parameters are evaluated:
# 0: radar-to-target distance beeline
# 1: projection of radar-to-target distance beeline on ground level
# 2: swath range
# 3: antenna footprint center
# 4,5,6,7: antenna footprint semiaxis (upper, lower, right, left)

tiltAngle_index = 0
direction_index = 0

for tiltAngle in tiltAngles:
    grid_parameters[tiltAngle_index, :, 0] = radarTargetDistance_projection = antennaHeight / np.sin(tiltAngle)
    grid_parameters[tiltAngle_index, :, 1] = antennaHeight / np.tan(tiltAngle)
    grid_parameters[tiltAngle_index, :, 2] = radarTargetDistance_projection * np.tan(MAX_SCAN_ANGLE) * 2
    print(f"Tilt Angle: {np.rad2deg(tiltAngle):.1f} degree")
    print(f"\tRadar-to-target distance: {radarTargetDistance_projection:.1f} m")
    print(f"\tRadar-to-target distance (projection): {grid_parameters[tiltAngle_index, 0, 1]:.1f} m")
    print(f"\tSwath range: {grid_parameters[tiltAngle_index, 0, 2]:.1f} m")
    for direction in directions:
        grid_parameters[tiltAngle_index, direction_index, 3] = radarTargetDistance_projection * np.tan(direction)
        grid_parameters[tiltAngle_index, direction_index, 4] = np.abs(antennaHeight * (1/np.tan(tiltAngle) - 1/np.tan(tiltAngle-BEAMWIDTH_ELEVATION/2)))
        grid_parameters[tiltAngle_index, direction_index, 5] = np.abs(antennaHeight * (1/np.tan(tiltAngle) - 1/np.tan(tiltAngle+BEAMWIDTH_ELEVATION/2)))
        grid_parameters[tiltAngle_index, direction_index, 6] = np.abs(radarTargetDistance_projection * np.tan(direction) - radarTargetDistance_projection * np.tan(direction+BEAMWIDTH_AZIMUTH/2))
        grid_parameters[tiltAngle_index, direction_index, 7] = np.abs(radarTargetDistance_projection * np.tan(direction) - radarTargetDistance_projection * np.tan(direction-BEAMWIDTH_AZIMUTH/2))
        plt.errorbar(grid_parameters[tiltAngle_index,direction_index,3], grid_parameters[tiltAngle_index,direction_index,1], xerr=[[grid_parameters[tiltAngle_index,direction_index,7]], [grid_parameters[tiltAngle_index,direction_index,6]]], yerr=[[grid_parameters[tiltAngle_index,direction_index,5]],[grid_parameters[tiltAngle_index,direction_index,4]]], fmt="o", color='r', markersize=7, capsize=10)
        direction_index += 1
    direction_index = 0
    tiltAngle_index += 1
plt.axis('equal')
plt.grid(visible=True, which="both", linewidth=0.3)
plt.show()

