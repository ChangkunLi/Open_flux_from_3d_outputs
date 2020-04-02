import numpy as np
import tecplot as tp
from numpy import abs, pi, cos, sin

RMercury = 2440000

# connect to Tecplot
tp.session.connect()

# read in data
tp.new_layout()
dataset = tp.data.load_tecplot('/Volumes/DATA_disk/Frontera_run/SW_test/RUN_46/3d__var_7_t00000118_n00692607.plt')
frame = tp.active_frame()

# creating spherical zone
shape = (32,32) # used for debugging
# shape = (128,128)

r = 1.05
phi = np.linspace(0, pi, shape[0])
theta = np.linspace(0, 2*pi, shape[1])

ttheta, pphi = np.meshgrid(theta, phi, indexing='ij')

xx = r * sin(pphi) * cos(ttheta)
yy = r * sin(pphi) * sin(ttheta)
zz = r * cos(pphi)

sphere_zone = dataset.add_ordered_zone('R = {}'.format(r),shape)

sphere_zone.values('X*')[:] = xx.ravel()
sphere_zone.values('Y*')[:] = yy.ravel()
sphere_zone.values('Z*')[:] = zz.ravel()

X_seed = sphere_zone.values('X*')
Y_seed = sphere_zone.values('Y*')
Z_seed = sphere_zone.values('Z*')

# interpolate magnetic field to spherical zone
tp.data.operate.interpolate_linear(sphere_zone, dataset.zone('3D*'), dataset.variable('B_x*'))
tp.data.operate.interpolate_linear(sphere_zone, dataset.zone('3D*'), dataset.variable('B_y*'))
tp.data.operate.interpolate_linear(sphere_zone, dataset.zone('3D*'), dataset.variable('B_z*'))

# set up vectors and background contour
plot = frame.plot()
plot.vector.u_variable = dataset.variable('B_x*')
plot.vector.v_variable = dataset.variable('B_y*')
plot.vector.w_variable = dataset.variable('B_z*')
plot.show_streamtraces = False

streamtraces = plot.streamtraces

streamtraces.step_size = .25
streamtraces.max_steps = 10000

# add status variable
dataset.add_variable('Status')

# initialize polar cap area and open flux
Area = {}
Flux = {}
Area['North'] = 0.0
Area['South'] = 0.0
Flux['North'] = 0.0
Flux['South'] = 0.0

# Pre-allocate theta value for each node
Phi = pphi.ravel()

# suspend the session
tp.session.suspend_enter()

# loop over nodes on 2D sphere_zone
for i in range(len(X_seed)):
    if np.mod(i+1,1000) == 0:
        print('Iteration {}'.format(i))
    streamtraces.add([X_seed[i], Y_seed[i], Z_seed[i]],tp.constant.Streamtrace.VolumeLine)
    slice_zone = streamtraces.extract()
    # slice_zone is a one-element generator, use a for loop to get its value
    for streamline in slice_zone:
        pass
    X = streamline.values('X*')
    Y = streamline.values('Y*')
    Z = streamline.values('Z*')
    
    X_max = np.fmax(np.abs(X.min()),np.abs(X.max()))
    Y_max = np.fmax(np.abs(Y.min()),np.abs(Y.max()))
    Z_max = np.fmax(np.abs(Z.min()),np.abs(Z.max()))
    if X_max < 3.0 and Y_max < 3.0 and Z_max < 3.0:
        status = 3       # close field line
    elif Z_seed[i] > 0.0:
        status = 2
    else:
        status = 1
    
    sphere_zone.values('Status')[i] = status

    # clear slice_zone and streamtraces
    dataset.delete_zones(dataset.zone('Streamtrace'))
    streamtraces.delete_all()

    # calculating open flux and area
    if status == 2:
        Area['North'] += r*r*sin(Phi[i])*(pi/shape[0])*(2*pi/shape[1])
        Bfield = np.array([sphere_zone.values('B_x*')[i], sphere_zone.values('B_y*')[i], sphere_zone.values('B_z*')[i]])
        normal = np.array([X_seed[i], Y_seed[i], Z_seed[i]])
        Flux['North'] += np.abs(np.dot(Bfield, normal))*r*sin(Phi[i])*(pi/shape[0])*(2*pi/shape[1])
    elif status == 1:
        Area['South'] += r*r*sin(Phi[i])*(pi/shape[0])*(2*pi/shape[1])
        Bfield = np.array([sphere_zone.values('B_x*')[i], sphere_zone.values('B_y*')[i], sphere_zone.values('B_z*')[i]])
        normal = np.array([X_seed[i], Y_seed[i], Z_seed[i]])
        Flux['South'] += np.abs(np.dot(Bfield, normal))*r*sin(Phi[i])*(pi/shape[0])*(2*pi/shape[1])

tp.session.suspend_exit()

# disconnect from Tecplot
tp.session.disconnect()
            
# units conversion
Area['North'] *= RMercury*RMercury
Area['South'] *= RMercury*RMercury
Flux['North'] *= RMercury*RMercury*1e-9
Flux['South'] *= RMercury*RMercury*1e-9

# output results to stdout or another function
print('Northern polar cap area is {} m*m'.format(Area['North']))
print('Southern polar cap area is {} m*m'.format(Area['South']))
print('Northern open flux is {} Wb'.format(Flux['North']))
print('Southern open flux is {} Wb'.format(Flux['South']))






