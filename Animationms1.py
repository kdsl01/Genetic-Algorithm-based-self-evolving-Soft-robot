import open3d as o3d  # http://www.open3d.org/
import numpy as np
import time
import csv
import numba as nb

Gravity = 9.81
dt = 0.0001
xtT = 0
u_s = 0.61
u_k = 0.47
damping = 0.9999


# Mass class structure,
# contain wiehgt position velocity and acceleration of the mass
class mass:
    def __init__(self, weights, positions, velocity, acceleration):
        self.m = weights
        self.p = positions
        self.v = velocity
        self.a = acceleration


class spring:
    def __init__(self, SpringC, Length, m1, m2):
        self.k = SpringC
        self.L0 = Length
        self.L0x = Length
        self.m1 = m1
        self.m2 = m2


def setView(ctr, camera_pos=(1, 1, 1), lookat=(0, 0, 0), up=(0, 0, 1)):
    """
    set the view given a view control handel ctr
    """
    ctr.set_constant_z_far(100)  # camera z far clip plane
    ctr.set_constant_z_near(0.1)  # camera z near clip plane


def massinit(weight, sp):
    peak = [0] * 8
    point = [0] * 8
    peak[0] = sp;
    peak[1] = [sp[0] + 0.1, sp[1], sp[2]]
    peak[2] = [sp[0] + 0.1, sp[1] + 0.1, sp[2]]
    peak[3] = [sp[0] + 0.1, sp[1] + 0.1, sp[2] + 0.1]
    peak[4] = [sp[0], sp[1] + 0.1, sp[2]]
    peak[5] = [sp[0], sp[1] + 0.1, sp[2] + 0.1]
    peak[6] = [sp[0], sp[1], sp[2] + 0.1]
    peak[7] = [sp[0] + 0.1, sp[1], sp[2] + 0.1]
    for i in range(8):
        point[i] = mass(weight, peak[i], [0, 0, 0], [0, 0, 0])
    return point


def massinitv3(weight, sp, spc):
    peak = [0] * 8
    point = [0] * 8
    peak[0] = sp
    for icct in range(1, 8):
        peak[icct] = [sp[0] + spc[icct - 1][0],sp[1] + spc[icct - 1][1],sp[2] + spc[icct - 1][2] ]
    for i in range(8):
        point[i] = mass(weight, peak[i], [0, 0, 0], [0, 0, 0])
    return point


def springinitv3(peak):
    lines = [0] * 28
    kcount = 0
    for i in range(8):
        for j in range(i + 1, 8):
            spc = 20000
            tp1 = i
            tp2 = j
            lengtht = np.sqrt(
                np.power(peak[tp1].p[0] - peak[tp2].p[0], 2) + np.power(peak[tp1].p[1] - peak[tp2].p[1], 2) + np.power(
                    peak[tp1].p[2] - peak[tp2].p[2], 2))
            lines[kcount] = spring(spc, lengtht, tp1, tp2)
            kcount +=1
    return lines


def springinit(peak):
    lines = [0] * 28
    # linepos = [[0, 1], [0, 2], [0, 3], [0, 4],
    #          [0, 5], [0, 6], [0, 7], [1, 2],
    #          [1, 3], [1, 4], [1, 6], [1, 7],
    #          [2, 3], [2, 4], [2, 5], [2, 6],
    #          [2, 7], [3, 4], [3, 5], [3, 6],
    #          [3, 7], [4, 5], [4, 6], [4, 7],
    #          [5, 6], [5, 7], [6, 7], [1, 5]]
    linepos = [[0, 1], [0, 6], [1, 7], [2, 3],
               [2, 4], [3, 5], [4, 5], [6, 7],
               [0, 4], [1, 2], [3, 7], [5, 6],
               [0, 2], [0, 5], [0, 7], [1, 3],
               [1, 4], [1, 6], [2, 5], [2, 7],
               [3, 4], [3, 6], [4, 6], [5, 7],
               [0, 3], [1, 5], [2, 6], [4, 7]]
    for j in range(28):
        spc = 20000
        tp1 = linepos[j][0]
        tp2 = linepos[j][1]
        lengtht = np.sqrt(
            np.power(peak[tp1].p[0] - peak[tp2].p[0], 2) + np.power(peak[tp1].p[1] - peak[tp2].p[1], 2) + np.power(
                peak[tp1].p[2] - peak[tp2].p[2], 2))
        lines[j] = spring(spc, lengtht, tp1, tp2)
    return lines


def customVisualization(geomtry_list):
    """
    helper function to create a visualization given a list of o3d geometry
    """
    vis = o3d.visualization.Visualizer()
    vis.create_window()
    for g in geomtry_list:
        vis.add_geometry(g)
    ctr = vis.get_view_control()
    setView(ctr)
    vis.run()
    vis.destroy_window()  # close the window when finished


def createPlane(r=0.5, dr=0.1):
    """
    return a plane located at (0,0,0),and with plane normal = (0,0,1)
    r: radius of the plane
    dr:discretization radius of the grid
    """
    bounds = np.array([[-r, -r, 0], [r, r, 0]]) / dr
    bounds = np.stack((np.floor(bounds[0]), np.ceil(bounds[1]))) * dr
    nx, ny, nz = np.ceil((bounds[1] - bounds[0]) / dr).astype(int)
    #     print(nx,ny)
    xyz = np.reshape([[[[i, j, 0], [i + 1, j, 0], [i, j + 1, 0],
                        [i, j + 1, 0], [i + 1, j, 0], [i + 1, j + 1, 0]] for i in range(nx - 1)] for j in
                      range(ny - 1)], (-1, 3))
    xyz = (xyz - ((nx - 1) / 2, (ny - 1) / 2, 0)) * dr
    #     xyz, bounds, (nx, ny, nz) = create_grid(bounds, dr)
    #     print(nx, ny, nz)
    triangles = np.arange(xyz.shape[0]).reshape((-1, 3))
    plane = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(
        xyz), o3d.utility.Vector3iVector(triangles))
    # assign checkerboard color pattern
    c0 = (0.323, 0.78, 0.321)  # first color
    c1 = (0.863, 0.62, 0.343)  # second color
    colors = np.reshape([[np.tile(c0 if (i + j) % 2 else c1, (6, 1)) for i in range(nx - 1)] for j in range(ny - 1)],
                        (-1, 3))
    plane.vertex_colors = o3d.utility.Vector3dVector(colors)
    plane.compute_triangle_normals()
    return plane


def add_mass(peak):
    points = [0] * 8
    for i in range(8):
        points[i] = peak[i].p
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(points)
    return pcd


def create_block(sp=[0, 0, 0]):
    peak = [0] * 8
    peak[0] = sp;  # 0,0,0
    peak[1] = [sp[0] + 0.1, sp[1], sp[2]]  # 0.1,0,0
    peak[2] = [sp[0] + 0.1, sp[1] + 0.1, sp[2]]  # 0.1,0.1,0
    peak[3] = [sp[0] + 0.1, sp[1] + 0.1, sp[2] + 0.1]  # 0.1,0.1,0.1
    peak[4] = [sp[0], sp[1] + 0.1, sp[2]]  # 0,0.1,0
    peak[5] = [sp[0], sp[1] + 0.1, sp[2] + 0.1]  # 0,0.1,0.1
    peak[6] = [sp[0], sp[1], sp[2] + 0.1]  # 0,0,0.1
    peak[7] = [sp[0] + 0.1, sp[1], sp[2] + 0.1]  # 0.1,0,0.1
    return peak


def draw_surface(peak):
    points = [0] * 8
    for i in range(8):
        points[i] = peak[i].p
    mesh = o3d.geometry.TriangleMesh()
    mesh.vertices = o3d.utility.Vector3dVector(points)
    np_triangles = np.array([[0, 1, 6],
                             [7, 6, 1],
                             [1, 2, 7],
                             [3, 7, 2],
                             [5, 6, 3],
                             [7, 3, 6],
                             [4, 0, 5],
                             [6, 5, 0],
                             [2, 4, 3],
                             [5, 3, 4],
                             [1, 0, 2],
                             [4, 2, 0]]).astype(np.int32)
    mesh.triangles = o3d.utility.Vector3iVector(np_triangles)
    return mesh


def drawcubev2(lines, peak, so):
    points = [0] * 8
    linesx = np.zeros((len(so), 2))
    for i in range(8):
        points[i] = peak[i].p
    for j in range(len(so)):
        linesx[j][0] = lines[int(so[j])].m1
        linesx[j][1] = lines[int(so[j])].m2

    color = [[0, 0, 0] for i in range(len(lines))]
    line_set = o3d.geometry.LineSet()
    line_set.points = o3d.utility.Vector3dVector(points)
    line_set.lines = o3d.utility.Vector2iVector(linesx)
    line_set.colors = o3d.utility.Vector3dVector(color)
    return line_set


def updatecube(line, peak, cT: float):
    fx = [0] * 8
    fy = [0] * 8
    fz = [0] * 8
    GPE = [0] * 8
    KE = [0] * 8
    SPE = [0] * 28
    GRE = [0] * 8
    V = [0] * 8
    for i in range(8):
        for j in range(28):
            if line[j].m1 == i or line[j].m2 == i:
                if line[j].m1 == i:
                    self = line[j].m1
                    opp = line[j].m2
                elif line[j].m2 == i:
                    self = line[j].m2
                    opp = line[j].m1
                X = np.power(peak[opp].p[0] - peak[self].p[0], 2)
                Y = np.power(peak[opp].p[1] - peak[self].p[1], 2)
                Z = np.power(peak[opp].p[2] - peak[self].p[2], 2)
                L = np.sqrt(X + Y + Z)
                F = line[j].k * (L - line[j].L0)
                fx[i] = fx[i] - F * ((peak[self].p[0] - peak[opp].p[0]) / L)
                fy[i] = fy[i] - F * ((peak[self].p[1] - peak[opp].p[1]) / L)
                fz[i] = fz[i] - F * ((peak[self].p[2] - peak[opp].p[2]) / L)
        fz[i] = fz[i] - peak[i].m * Gravity
        if peak[i].p[2] <= 0:
            fz[i] = fz[i] + (0 - peak[i].p[2]) * 100000
            GRE[i] = 0.5 * 100000 * np.power(peak[i].p[2], 2)
    for count in range(8):
        peak[count].a[0] = fx[count] / peak[count].m
        peak[count].a[1] = fy[count] / peak[count].m
        peak[count].a[2] = fz[count] / peak[count].m
        peak[count].v[0] = peak[count].v[0] + peak[count].a[0] * dt
        peak[count].v[1] = peak[count].v[1] + peak[count].a[1] * dt
        peak[count].v[2] = peak[count].v[2] + peak[count].a[2] * dt
        peak[count].v[0] = damping * peak[count].v[0]
        peak[count].v[1] = damping * peak[count].v[1]
        peak[count].v[2] = damping * peak[count].v[2]
        peak[count].p[0] = peak[count].p[0] + peak[count].v[0] * dt
        peak[count].p[1] = peak[count].p[1] + peak[count].v[1] * dt
        peak[count].p[2] = peak[count].p[2] + peak[count].v[2] * dt
        GPE[count] = peak[count].p[2] * peak[count].m * Gravity
        V[count] = np.sqrt(
            np.power(peak[count].v[0], 2) + np.power(peak[count].v[1], 2) + np.power(peak[count].v[2], 2))
        KE[count] = 0.5 * peak[count].m * np.power(V[count], 2)
    for countx in range(28):
        xx = np.power(peak[line[countx].m1].p[0] - peak[line[countx].m2].p[0], 2)
        yy = np.power(peak[line[countx].m1].p[1] - peak[line[countx].m2].p[1], 2)
        zz = np.power(peak[line[countx].m1].p[2] - peak[line[countx].m2].p[2], 2)
        LL = np.sqrt(xx + yy + zz)
        SPE[countx] = 0.5 * line[countx].k * np.power(LL - line[countx].L0, 2)
    cT += dt
    return cT, np.sum(GPE), np.sum(KE), np.sum(SPE), np.sum(GRE)


def updatecube_friction(line, peak, cT, so):
    fx = [0] * 8
    fy = [0] * 8
    fz = [0] * 8
    groundflag = [0] * 8
    for flagcount in range(8):
        groundflag[flagcount] = False
    for i in range(8):
        for j in range(len(so)):
            if line[int(so[j])].m1 == i or line[int(so[j])].m2 == i:
                if line[int(so[j])].m1 == i:
                    self = line[int(so[j])].m1
                    opp = line[int(so[j])].m2
                elif line[int(so[j])].m2 == i:
                    self = line[int(so[j])].m2
                    opp = line[int(so[j])].m1
                X = np.power(peak[opp].p[0] - peak[self].p[0], 2)
                Y = np.power(peak[opp].p[1] - peak[self].p[1], 2)
                Z = np.power(peak[opp].p[2] - peak[self].p[2], 2)
                L = np.sqrt(X + Y + Z)
                F = line[int(so[j])].k * (L - line[int(so[j])].L0)
                fx[i] = fx[i] - F * ((peak[self].p[0] - peak[opp].p[0]) / L)
                fy[i] = fy[i] - F * ((peak[self].p[1] - peak[opp].p[1]) / L)
                fz[i] = fz[i] - F * ((peak[self].p[2] - peak[opp].p[2]) / L)
        fz[i] = fz[i] - peak[i].m * Gravity
        if peak[i].p[2] <= 0:
            Fh = np.sqrt(np.power(fx[i], 2) + np.power(fy[i], 2))
            if Fh < abs(fz[i] * u_s):
                groundflag[i] = True
            else:
                fx[i] = fx[i] - (fx[i] / Fh) * abs(fz[i] * u_k)
                fy[i] = fy[i] - (fy[i] / Fh) * abs(fz[i] * u_k)
            fz[i] = fz[i] + (0 - peak[i].p[2]) * 100000
    for count in range(8):
        peak[count].a[0] = fx[count] / peak[count].m
        peak[count].a[1] = fy[count] / peak[count].m
        peak[count].a[2] = fz[count] / peak[count].m
        peak[count].v[0] = peak[count].v[0] + peak[count].a[0] * dt
        peak[count].v[1] = peak[count].v[1] + peak[count].a[1] * dt
        peak[count].v[2] = peak[count].v[2] + peak[count].a[2] * dt
        peak[count].v[0] = damping * peak[count].v[0]
        peak[count].v[1] = damping * peak[count].v[1]
        peak[count].v[2] = damping * peak[count].v[2]
        if groundflag[count] == True:
            peak[count].v[0] = 0
            peak[count].v[1] = 0
            groundflag[count] = False
        peak[count].p[0] = peak[count].p[0] + peak[count].v[0] * dt
        peak[count].p[1] = peak[count].p[1] + peak[count].v[1] * dt
        peak[count].p[2] = peak[count].p[2] + peak[count].v[2] * dt
    cT += dt
    return cT


def breath(line, T, b, c):
    for count in range(0, 12):
        line[count].L0 = line[count].L0x + np.sin(23 * T + c[count]) / b[count]
    for count in range(12, 24):
        line[count].L0 = line[count].L0x + np.sin(23 * T + c[count]) / b[count]
    for count in range(24, 28):
        line[count].L0 = line[count].L0x + np.sin(23 * T + c[count]) / b[count]
    # print(T)
    # line[1].L0 = 0.1 + np.sin(23 * T) / 15
    # line[5].L0 = 0.1 + np.sin(23 * T) / 15
    # line[11].L0 = 0.1 + np.sin(23 * T) / 15
    # line[12].L0 = 0.1 + np.sin(23 * T) / 15
    # line[9].L0 = 0.1 + np.sin(23 * T) / 15
    # line[19].L0 = 0.1 + np.sin(23 * T) / 15
    # line[21].L0 = 0.1 + np.sin(23 * T) / 15
    # line[25].L0 = 0.1 + np.sin(23 * T) / 15
    # line[18].L0 = 0.1 + np.sin(23 * T) / 15
    # line[26].L0 = 0.1 + np.sin(23 * T) / 15
    #
    # line[5].L0 = 0.1 + np.sin(23 * T) / 15
    # line[11].L0 = 0.1 + np.sin(23 * T) / 15
    # line[12].L0 = 0.1 + np.sin(23 * T) / 15
    # line[21].L0 = 0.1 + np.sin(23 * T) / 15


def random_g():
    b = [0] * 28
    c = [0] * 28
    for i in range(0, 12):
        b[i] = np.random.uniform(30, 100)
        c[i] = np.random.uniform(0, 19)
    for j in range(12, 24):
        b[j] = np.random.uniform(45, 100)
        c[j] = np.random.uniform(0, 19)
    for k in range(24, 28):
        b[k] = np.random.uniform(44, 100)
        c[k] = np.random.uniform(0, 19)
    return b, c


def main():
    csv_file = open('GA_Robot_ms2v2.csv', 'r')
    a = []
    b = []
    c = []
    d = []
    sk = []

    # Read off and discard first line, to skip headers
    # Split columns while reading
    for dataA, dataB, dataC, dataD, dataE in csv.reader(csv_file, delimiter=','):
        # Append each variable to a separate list
        a.append(float(dataA))
        d.append(float(dataB))
        b.append(float(dataC))
        c.append(float(dataD))
        sk.append(int(dataE))
    # csv_file2 = open('GA_Robot_ms1_position.csv', 'r')
    # point_position = [];
    # for pa, pb, pc in csv.reader(csv_file2, delimiter=','):
    #     point_position.append([float(pa),float(pb),float(pc)])
    csv_file3 = open('springorder_ms2v2.csv','r')
    springorder = []
    for so in csv.reader(csv_file3, delimiter=','):
        springorder.append(so[0])
    T = 0
    ended = False
    pts = [0] * 8
    ptsx = [0] * 8

    def signalEnd(vis):
        nonlocal ended
        ended = True

    vis = o3d.visualization.VisualizerWithKeyCallback()
    plane = createPlane()
    peak = massinit(1, [0.1, 0.1, 0.1])
    coord_frame = o3d.geometry.TriangleMesh.create_coordinate_frame(size=0.1, origin=[0, 0, 0])
    c1 = []
    c2 = []
    c3 = []
    for xtc in peak:
        c1.append(xtc.p[0])
        c2.append(xtc.p[1])
        c3.append(xtc.p[2])
    ctm = [np.mean(c1), np.mean(c2), np.mean(c3)]
    print(ctm)
    print(peak[5].p)
    edges = springinitv3(peak)
    sf = draw_surface(peak)
    for i in range(len(edges)):
        edges[i].k = sk[i];
    lsd = drawcubev2(edges, peak,springorder)
    vis.register_key_callback(81, signalEnd)  # key Q
    vis.register_key_callback(256, signalEnd)  # key escape
    vis.create_window()
    vis.add_geometry(lsd)
    vis.add_geometry(plane)
    vis.add_geometry(coord_frame)
    ctr = vis.get_view_control()
    setView(ctr)
    t_prev = time.time()
    t_btv = time.time()

    xtxxx = 0
    xtxxxx = 0
    bflag = True
    # b, c = random_g()
    print(len(b))
    t_start = time.time()
    while (not ended):
        t = time.time()
        xtxxx = updatecube_friction(edges, peak, xtxxx,springorder)
        # breath(edges, xtxxx, b, c)
        if t - t_prev > 1. / 120.:
            t_prev = t
            for count in range(8):
                pts[count] = peak[count].p
            lsd.points = o3d.utility.Vector3dVector(pts)
            vis.update_geometry(lsd)
            vis.poll_events()
            vis.update_renderer()
        T = T + 1
    runtime = time.time() - t_start
    c1f = []
    c2f = []
    c3f = []
    for xtcf in peak:
        c1f.append(xtcf.p[0])
        c2f.append(xtcf.p[1])
        c3f.append(xtcf.p[2])
    ctmf = [np.mean(c1f), np.mean(c2f), np.mean(c3f)]
    dist = np.sqrt(np.power(ctmf[0] - ctm[0], 2) + np.power(ctmf[1] - ctm[1], 2))
    speed = dist / runtime
    print(runtime)
    print(ctmf)
    print(dist)
    print(speed)
    print(peak[5].p)
    print(dist / xtxxx)

    # customVisualization([lsd,spheres,plane])


main()

print("pass")