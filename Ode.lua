local ffi = require("ffi")

-- https://ode.org/wiki/index.php/Main_Page
-- https://bitbucket.org/odedevs/ode/src/master/include/ode/

-- odeconfig.h
ffi.cdef[[
]]

-- common.h
ffi.cdef[[
    typedef double dReal; //57

    typedef enum {
        dSA__MIN,
    
        dSA_X = dSA__MIN,
        dSA_Y,
        dSA_Z,
    
        dSA__MAX,
    } dSpaceAxis; //93

    typedef enum {
        dMD__MIN,
    
        dMD_LINEAR = dMD__MIN,
        dMD_ANGULAR,
    
        dMD__MAX,
    } dMotionDynamics; //110

    typedef enum {
        dDA__MIN,
    
        dDA__L_MIN = dDA__MIN + dMD_LINEAR * dSA__MAX,
    
        dDA_LX = dDA__L_MIN + dSA_X,
        dDA_LY = dDA__L_MIN + dSA_Y,
        dDA_LZ = dDA__L_MIN + dSA_Z,
    
        dDA__L_MAX = dDA__L_MIN + dSA__MAX,
    
        dDA__A_MIN = dDA__MIN + dMD_ANGULAR * dSA__MAX,
    
        dDA_AX = dDA__A_MIN + dSA_X,
        dDA_AY = dDA__A_MIN + dSA_Y,
        dDA_AZ = dDA__A_MIN + dSA_Z,
    
        dDA__A_MAX = dDA__A_MIN + dSA__MAX,
    
        dDA__MAX = dDA__MIN + dMD__MAX * dSA__MAX,
    } dDynamicsAxis; //112

    typedef enum {
        dV3E__MIN,
    
        dV3E__AXES_MIN = dV3E__MIN,
    
        dV3E_X = dV3E__AXES_MIN + dSA_X,
        dV3E_Y = dV3E__AXES_MIN + dSA_Y,
        dV3E_Z = dV3E__AXES_MIN + dSA_Z,
    
        dV3E__AXES_MAX = dV3E__AXES_MIN + dSA__MAX,
    
        dV3E_PAD = dV3E__AXES_MAX,
    
        dV3E__MAX,
    
        dV3E__AXES_COUNT = dV3E__AXES_MAX - dV3E__AXES_MIN,
    } dVec3Element; //134

    typedef enum {
        dV4E__MIN,
    
        dV4E_X = dV4E__MIN + dSA_X,
        dV4E_Y = dV4E__MIN + dSA_Y,
        dV4E_Z = dV4E__MIN + dSA_Z,
        dV4E_O = dV4E__MIN + dSA__MAX,
    
        dV4E__MAX,
    } dVec4Element; //152

    typedef enum {
        dM3E__MIN,
    
        dM3E__X_MIN = dM3E__MIN + dSA_X * dV3E__MAX,
        
        dM3E__X_AXES_MIN = dM3E__X_MIN + dV3E__AXES_MIN,
    
        dM3E_XX = dM3E__X_MIN + dV3E_X,
        dM3E_XY = dM3E__X_MIN + dV3E_Y,
        dM3E_XZ = dM3E__X_MIN + dV3E_Z,
    
        dM3E__X_AXES_MAX = dM3E__X_MIN + dV3E__AXES_MAX,
    
        dM3E_XPAD = dM3E__X_MIN + dV3E_PAD,
    
        dM3E__X_MAX = dM3E__X_MIN + dV3E__MAX,
    
        dM3E__Y_MIN = dM3E__MIN + dSA_Y * dV3E__MAX,
    
        dM3E__Y_AXES_MIN = dM3E__Y_MIN + dV3E__AXES_MIN,
    
        dM3E_YX = dM3E__Y_MIN + dV3E_X,
        dM3E_YY = dM3E__Y_MIN + dV3E_Y,
        dM3E_YZ = dM3E__Y_MIN + dV3E_Z,
    
        dM3E__Y_AXES_MAX = dM3E__Y_MIN + dV3E__AXES_MAX,
    
        dM3E_YPAD = dM3E__Y_MIN + dV3E_PAD,
    
        dM3E__Y_MAX = dM3E__Y_MIN + dV3E__MAX,
    
        dM3E__Z_MIN = dM3E__MIN + dSA_Z * dV3E__MAX,
    
        dM3E__Z_AXES_MIN = dM3E__Z_MIN + dV3E__AXES_MIN,
    
        dM3E_ZX = dM3E__Z_MIN + dV3E_X,
        dM3E_ZY = dM3E__Z_MIN + dV3E_Y,
        dM3E_ZZ = dM3E__Z_MIN + dV3E_Z,
    
        dM3E__Z_AXES_MAX = dM3E__Z_MIN + dV3E__AXES_MAX,
    
        dM3E_ZPAD = dM3E__Z_MIN + dV3E_PAD,
    
        dM3E__Z_MAX = dM3E__Z_MIN + dV3E__MAX,
    
        dM3E__MAX = dM3E__MIN + dSA__MAX * dV3E__MAX,
    } dMat3Element; //209

    typedef enum {
        dM4E__MIN,
    
        dM4E__X_MIN = dM4E__MIN + dV4E_X * dV4E__MAX,
    
        dM4E_XX = dM4E__X_MIN + dV4E_X,
        dM4E_XY = dM4E__X_MIN + dV4E_Y,
        dM4E_XZ = dM4E__X_MIN + dV4E_Z,
        dM4E_XO = dM4E__X_MIN + dV4E_O,
    
        dM4E__X_MAX = dM4E__X_MIN + dV4E__MAX,
    
        dM4E__Y_MIN = dM4E__MIN + dV4E_Y * dV4E__MAX,
    
        dM4E_YX = dM4E__Y_MIN + dV4E_X,
        dM4E_YY = dM4E__Y_MIN + dV4E_Y,
        dM4E_YZ = dM4E__Y_MIN + dV4E_Z,
        dM4E_YO = dM4E__Y_MIN + dV4E_O,
    
        dM4E__Y_MAX = dM4E__Y_MIN + dV4E__MAX,
    
        dM4E__Z_MIN = dM4E__MIN + dV4E_Z * dV4E__MAX,
    
        dM4E_ZX = dM4E__Z_MIN + dV4E_X,
        dM4E_ZY = dM4E__Z_MIN + dV4E_Y,
        dM4E_ZZ = dM4E__Z_MIN + dV4E_Z,
        dM4E_ZO = dM4E__Z_MIN + dV4E_O,
    
        dM4E__Z_MAX = dM4E__Z_MIN + dV4E__MAX,
    
        dM4E__O_MIN = dM4E__MIN + dV4E_O * dV4E__MAX,
    
        dM4E_OX = dM4E__O_MIN + dV4E_X,
        dM4E_OY = dM4E__O_MIN + dV4E_Y,
        dM4E_OZ = dM4E__O_MIN + dV4E_Z,
        dM4E_OO = dM4E__O_MIN + dV4E_O,
    
        dM4E__O_MAX = dM4E__O_MIN + dV4E__MAX,
    
        dM4E__MAX = dM4E__MIN + dV4E__MAX * dV4E__MAX,
    } dMat4Element; //211

    typedef enum {
        dQUE__MIN,
    
        dQUE_R = dQUE__MIN,
    
        dQUE__AXIS_MIN,
    
        dQUE_I = dQUE__AXIS_MIN + dSA_X,
        dQUE_J = dQUE__AXIS_MIN + dSA_Y,
        dQUE_K = dQUE__AXIS_MIN + dSA_Z,
    
        dQUE__AXIS_MAX = dQUE__AXIS_MIN + dSA__MAX,
    
        dQUE__MAX = dQUE__AXIS_MAX,
    } dQuatElement; //253

    typedef dReal dVector3[dV3E__MAX]; //270
    typedef dReal dVector4[dV4E__MAX]; //271
    typedef dReal dMatrix3[dM3E__MAX]; //272
    typedef dReal dMatrix4[dM4E__MAX]; //273
    typedef dReal dMatrix6[(dMD__MAX * dV3E__MAX) * (dMD__MAX * dSA__MAX)]; //274
    typedef dReal dQuaternion[dQUE__MAX]; //275

    typedef struct dxWorld *dWorldID; //364
    typedef struct dxSpace *dSpaceID;
    typedef struct dxBody *dBodyID;
    typedef struct dxGeom *dGeomID;
    typedef struct dxJoint *dJointID;
    typedef struct dxJointGroup *dJointGroupID;

    enum {
        dLimotLoStop		= 0x0001,
        dLimotHiStop		= 0x0002,
        dLimotVel		= 0x0004,
        dLimotFMax		= 0x0008,
        dLimotFudgeFactor	= 0x0010,
        dLimotBounce		= 0x0020,
        dLimotSoft		= 0x0040
    }; //420

    enum {
        dAMotorUser = 0,
        dAMotorEuler = 1
    }; //500

    typedef struct dJointFeedback {
        dVector3 f1;		/* force applied to body 1 */
        dVector3 t1;		/* torque applied to body 1 */
        dVector3 f2;		/* force applied to body 2 */
        dVector3 t2;		/* torque applied to body 2 */
    } dJointFeedback; //516
]]

-- odeinit.h
ffi.cdef[[
    void dInitODE(void);
    void dCloseODE(void);
]]

-- contact.h
ffi.cdef[[
    enum {
        dContactMu2	  = 0x001,      /**< Use axis dependent friction */
        dContactAxisDep = 0x001,      /**< Same as above */
        dContactFDir1	  = 0x002,      /**< Use FDir for the first friction value */
        dContactBounce  = 0x004,      /**< Restore collision energy anti-parallel to the normal */
        dContactSoftERP = 0x008,      /**< Don't use global erp for penetration reduction */
        dContactSoftCFM = 0x010,      /**< Don't use global cfm for penetration constraint */
        dContactMotion1 = 0x020,      /**< Use a non-zero target velocity for the constraint */
        dContactMotion2 = 0x040, 
        dContactMotionN = 0x080, 
        dContactSlip1	  = 0x100,      /**< Force-dependent slip. */
        dContactSlip2	  = 0x200, 
        dContactRolling = 0x400,      /**< Rolling/Angular friction */
      
        dContactApprox0   = 0x0000,
        dContactApprox1_1 = 0x1000,
        dContactApprox1_2 = 0x2000,
        dContactApprox1_N = 0x4000,   /**< For rolling friction */
        dContactApprox1   = 0x7000
      }; //33

      typedef struct dSurfaceParameters {
        /* must always be defined */
        int   mode;
        dReal mu;
      
        /* only defined if the corresponding flag is set in mode */
        dReal mu2;
        dReal rho;                    /**< Rolling friction */
        dReal rho2;
        dReal rhoN;                   /**< Spinning friction */
        dReal bounce;                 /**< Coefficient of restitution */
        dReal bounce_vel;             /**< Bouncing threshold */
        dReal soft_erp;
        dReal soft_cfm;
        dReal motion1,motion2,motionN;
        dReal slip1,slip2;
      } dSurfaceParameters; //55

    typedef struct dContactGeom {
        dVector3 pos;          /*< contact position*/
        dVector3 normal;       /*< normal vector*/
        dReal depth;           /*< penetration depth*/
        dGeomID g1,g2;         /*< the colliding geoms*/
        int side1,side2;       /*< (to be documented)*/
    } dContactGeom; //88

    typedef struct dContact {
        dSurfaceParameters surface;
        dContactGeom geom;
        dVector3 fdir1;
      } dContact; //99
]]

-- mass.h
ffi.cdef[[
    struct dMass;
    typedef struct dMass dMass;

    void dMassSetSphere (dMass *, dReal density, dReal radius); //52
    void dMassSetCapsule (dMass *, dReal density, int direction, dReal radius, dReal length); //55
    void dMassSetBox (dMass *, dReal density, dReal lx, dReal ly, dReal lz); //65

    struct dMass {
        dReal mass;
        dVector3 c;
        dMatrix3 I;
    }; //88
]]

-- collision_space.h
ffi.cdef[[
    struct dContactGeom; //32
    typedef void (*dNearCallback) (void *data, dGeomID o1, dGeomID o2); //49 CHANGE!!!
    dSpaceID dSimpleSpaceCreate (dSpaceID space); //52
    dSpaceID dHashSpaceCreate (dSpaceID space); //53
    void dSpaceDestroy (dSpaceID); //70
    void dSpaceAdd (dSpaceID, dGeomID); //149
    void dSpaceRemove (dSpaceID, dGeomID); //150
]]

-- collision.h
ffi.cdef[[
    void dGeomDestroy (dGeomID geom); //65
    void dGeomSetBody (dGeomID geom, dBodyID body); //105
    dBodyID dGeomGetBody (dGeomID geom); //114
    void dGeomSetPosition (dGeomID geom, dReal x, dReal y, dReal z); //131
    void dGeomSetRotation (dGeomID geom, const dMatrix3 R); //146
    void dGeomSetQuaternion (dGeomID geom, const dQuaternion Q); //162
    const dReal * dGeomGetPosition (dGeomID geom); //181
    const dReal * dGeomGetRotation (dGeomID geom); //210
    int dCollide (dGeomID o1, dGeomID o2, int flags, dContactGeom *contact, int skip); //792
    void dSpaceCollide (dSpaceID space, void *data, dNearCallback callback); //822 CHANGE!!
    dGeomID dCreateSphere (dSpaceID space, dReal radius); //921
    dGeomID dCreateBox (dSpaceID space, dReal lx, dReal ly, dReal lz); //1001
    void dGeomBoxSetLengths (dGeomID box, dReal lx, dReal ly, dReal lz); //1015
    void dGeomBoxGetLengths (dGeomID box, dVector3 result); //1027
    dGeomID dCreatePlane (dSpaceID space, dReal a, dReal b, dReal c, dReal d); //1045
    void dGeomPlaneSetParams (dGeomID plane, dReal a, dReal b, dReal c, dReal d);
    void dGeomPlaneGetParams (dGeomID plane, dVector4 result);
    dReal dGeomPlanePointDepth (dGeomID plane, dReal x, dReal y, dReal z);
    dGeomID dCreateCapsule (dSpaceID space, dReal radius, dReal length); //1050
    void dGeomCapsuleSetParams (dGeomID ccylinder, dReal radius, dReal length); //1051
    void dGeomCapsuleGetParams (dGeomID ccylinder, dReal *radius, dReal *length); //1052
    dReal dGeomCapsulePointDepth (dGeomID ccylinder, dReal x, dReal y, dReal z); //1053
    dGeomID dCreateRay (dSpaceID space, dReal length); //1066
    void dGeomRaySetLength (dGeomID ray, dReal length); //1067
    dReal dGeomRayGetLength (dGeomID ray); //1068
    void dGeomRaySet (dGeomID ray, dReal px, dReal py, dReal pz, dReal dx, dReal dy, dReal dz); //1069
    void dGeomRayGet (dGeomID ray, dVector3 start, dVector3 dir); //1071
]]

-- objects.h
ffi.cdef[[
    // World
    dWorldID dWorldCreate(void); //52
    void dWorldDestroy(dWorldID world); //64
    void dWorldSetData (dWorldID world, void* data); //73
    void* dWorldGetData (dWorldID world); //82
    void dWorldSetGravity (dWorldID, dReal x, dReal y, dReal z); //93
    void dWorldGetGravity (dWorldID, dVector3 gravity); //100
    void dWorldSetERP (dWorldID, dReal erp); //110
    dReal dWorldGetERP (dWorldID); //117
    void dWorldSetCFM (dWorldID, dReal cfm); //127
    dReal dWorldGetCFM (dWorldID); //134
    
    int dWorldStep (dWorldID w, dReal stepsize); //384
    int dWorldQuickStep (dWorldID w, dReal stepsize); //427
    void dWorldSetLinearDamping (dWorldID w, dReal scale); //826
    void dWorldSetAngularDamping(dWorldID w, dReal scale); //840

    // Body
    dWorldID dBodyGetWorld (dBodyID); //1003
    dBodyID dBodyCreate (dWorldID); //1011
    void dBodyDestroy (dBodyID); //1021
    void dBodySetPosition   (dBodyID, dReal x, dReal y, dReal z); //1045
    const dReal * dBodyGetPosition (dBodyID); //1088
    const dReal * dBodyGetRotation (dBodyID); //1106
    const dReal * dBodyGetQuaternion (dBodyID); //1124
    void dBodySetMass (dBodyID, const dMass *mass); //1053
    void dBodyGetMass (dBodyID, dMass *mass); //1059
    void dBodySetQuaternion (dBodyID, const dQuaternion q); //1065
    void dBodySetLinearVel  (dBodyID, dReal x, dReal y, dReal z); //1071
    void dBodySetAngularVel (dBodyID, dReal x, dReal y, dReal z); //1077
    const dReal * dBodyGetLinearVel (dBodyID); //1141
    const dReal * dBodyGetAngularVel (dBodyID); //1147
    void dBodyAddForce (dBodyID, dReal fx, dReal fy, dReal fz); //1165
    void dBodyAddTorque (dBodyID, dReal fx, dReal fy, dReal fz); //1171
    const dReal * dBodyGetForce (dBodyID); //1219
    const dReal * dBodyGetTorque (dBodyID); //1230
    void dBodySetForce  (dBodyID b, dReal x, dReal y, dReal z); //1240
    void dBodySetTorque (dBodyID b, dReal x, dReal y, dReal z); //1250
    void dBodySetGravityMode (dBodyID b, int mode); //1449
    int dBodyGetGravityMode (dBodyID b); //1456
    dReal dBodyGetLinearDamping (dBodyID b); //1506
    void dBodySetLinearDamping(dBodyID b, dReal scale); //1516
    dReal dBodyGetAngularDamping (dBodyID b); //1524
    void dBodySetAngularDamping(dBodyID b, dReal scale); //1534
    void dBodySetDamping(dBodyID b, dReal linear_scale, dReal angular_scale); //1543
    dReal dBodyGetLinearDampingThreshold (dBodyID b); //1549
    void dBodySetLinearDampingThreshold(dBodyID b, dReal threshold); //1557
    dReal dBodyGetAngularDampingThreshold (dBodyID b); //1563
    void dBodySetAngularDampingThreshold(dBodyID b, dReal threshold); //1571
    dReal dBodyGetMaxAngularSpeed (dBodyID b); //1578
    void dBodySetMaxAngularSpeed(dBodyID b, dReal max_speed); //1588

    // Joint
    dJointID dJointCreateBall (dWorldID, dJointGroupID); //1693
    dJointID dJointCreateContact (dWorldID, dJointGroupID, const dContact *); //1717
    dJointID dJointCreateFixed (dWorldID, dJointGroupID); //1766
    dJointID dJointCreateAMotor (dWorldID, dJointGroupID); //1776
    dJointID dJointCreateLMotor (dWorldID, dJointGroupID); //1784
    void dJointDestroy (dJointID); //1827
    dJointGroupID dJointGroupCreate (int max_size); //1835
    void dJointGroupDestroy (dJointGroupID); //1843
    void dJointGroupEmpty (dJointGroupID); //1852
    void dJointAttach (dJointID, dBodyID body1, dBodyID body2); //1873
    void dJointSetFeedback (dJointID, dJointFeedback *);
    dJointFeedback *dJointGetFeedback (dJointID); //1956
    void dJointSetBallAnchor (dJointID, dReal x, dReal y, dReal z); //1965
    void dJointSetBallAnchor2 (dJointID, dReal x, dReal y, dReal z); //1971
    void dJointSetBallParam (dJointID, int parameter, dReal value); //1977
    void dJointSetAMotorNumAxes (dJointID, int num); //2463
    void dJointSetAMotorAxis (dJointID, int anum, int rel, dReal x, dReal y, dReal z); //2469
    void dJointSetAMotorAngle (dJointID, int anum, dReal angle); //2481
    void dJointSetAMotorParam (dJointID, int parameter, dReal value); //2487
    void dJointSetAMotorMode (dJointID, int mode); //2493
    void dJointSetFixed (dJointID); //2499
    void dJointSetFixedParam (dJointID, int parameter, dReal value); //2456
    void dJointAddAMotorTorques (dJointID, dReal torque1, dReal torque2, dReal torque3); //2503
    int dJointGetAMotorMode (dJointID); //3118
]]


local path = "/home/diminiminim/GameDev/VSLike/src/libode.so"
local Ode = ffi.load(path)

--[[
    local world = Ode.dWorldCreate()
    Ode.dWorldSetGravity(world, 0,-10,0)
    Ode.dWorldSetLinearDamping(world, 0.02)
    Ode.dWorldSetAngularDamping(world, 0.02)
    Ode.dWorldSetERP(world, 0.2)
    Ode.dWorldSetCFM(world, 1e-5)

    local space = Ode.dSimpleSpaceCreate(nil)

    local body = Ode.dBodyCreate(world)
    local geom = Ode.dCreateBox(space, 1,1,1)
    local mass = ffi.new("dMass")

    Ode.dMassSetBox (mass, 1, 3,3,3);
    Ode.dBodySetMass(body, mass);
    Ode.dGeomSetBody(geom, body)
    Ode.dBodySetPosition(body, set.x,set.y,set.z)

    Ode.dCreatePlane(space, 0,1,0,0)

    Ode.dBodyAddForce(body, dx, dy, dz)
    Ode.dBodyAddTorque(body, 0, dry, 0)

    local contactGroup = Ode.dJointGroupCreate(0)

    local function callback(data, o1, o2)
        local b1 = Ode.dGeomGetBody(o1)
        local b2 = Ode.dGeomGetBody(o2)
        local world = Ode.dBodyGetWorld(b1)
        local contact = ffi.new("dContact")

        
        if Ode.dCollide(o1, o2, 1, contact.geom, ffi.sizeof(contact)) then
            local joint = Ode.dJointCreateContact(world, contactGroup, contact)

            Ode.dJointAttach(joint, b1, b2)
        end
    end
        
    callback_C_func = ffi.cast("dNearCallback", callback)

    Ode.dSpaceCollide(e.odeSpace._, nil, callback_C_func)
    callback_C_func:free()

    Ode.dWorldStep(e.odeWorld._, 1/60)

    Ode.dJointGroupEmpty(contactGroup);

]]

return {Ode, ffi}