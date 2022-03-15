from tracemalloc import stop
import numpy as np

PIVOT_HEIGHT = 0.56 + 3.0 # m. Vertical distance between pivot and ground level.
ANTENNA_CENTER_POSITION = 0.12 # m. Distance between pivot and center of RX antenna.
MAX_SCAN_ANGLE = np.deg2rad(15) # rad
INITIAL_TILT_ANGLE = np.deg2rad(15) # rad
FINAL_TILT_ANGLE = np.deg2rad(45) # rad
TILT_POSITIONS = 3
SCANNING_DIRECTIONS = 3
BEAMWIDTH_AZIMUTH = np.deg2rad(6) # rad
BEAMWIDTH_ELEVATION = np.deg2rad(9) # rad


antennaHeight = PIVOT_HEIGHT + ANTENNA_CENTER_POSITION * np.cos(INITIAL_TILT_ANGLE)
tiltAngles = np.linspace(start=INITIAL_TILT_ANGLE, stop=FINAL_TILT_ANGLE, num=TILT_POSITIONS, endpoint=True)
directions = np.linspace(start=MAX_SCAN_ANGLE, stop=-MAX_SCAN_ANGLE, num=SCANNING_DIRECTIONS, endpoint=True)

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

for tiltAngle in tiltAngles:
    radarTargetDistance_beeline = antennaHeight / np.sin(tiltAngle)
    radarTargetDistance_projection = antennaHeight / np.tan(tiltAngle)
    swathRange = radarTargetDistance_projection * np.tan(MAX_SCAN_ANGLE) * 2
    for direction in directions:
        antennaFootprint_swathPosition = radarTargetDistance_projection * np.tan(direction)
        antennaFootprint_upperSemiaxis = np.abs(antennaHeight * (1/np.tan(tiltAngle) - 1/np.tan(tiltAngle-BEAMWIDTH_ELEVATION/2)))
        antennaFootprint_lowerSemiaxis = np.abs(antennaHeight * (1/np.tan(tiltAngle) - 1/np.tan(tiltAngle+BEAMWIDTH_ELEVATION/2)))
        antennaFootprint_rightSemiaxis = np.abs(radarTargetDistance_projection * np.tan(direction) - radarTargetDistance_projection * np.tan(direction+BEAMWIDTH_AZIMUTH/2))
        antennaFootprint_leftSemiaxis = np.abs(radarTargetDistance_projection * np.tan(direction) - radarTargetDistance_projection * np.tan(direction-BEAMWIDTH_AZIMUTH/2))