#include <cctk.h>
#include <cctk_Schedule.h>
#include <cctki_ScheduleBindings.h>
#include <cctki_GHExtensions.h>
#include <cctk_Parameters.h>
#include <iostream>
#include <assert.h>
#include <map>
#include <vector>

#define HERE { std::cout << "here " << __FILE__ << ":" << __LINE__ << std::endl; }

const char *where = "MyEvolve";

const char *ExtensionName = "MyDriver";

struct GF {
};

struct MyDriverData {
    std::map<int,int> time_levels;
};

int GroupStorageIncrease(const cGH *cctkGH, int n_groups, const int *groups,
                       const int *requested_tls, int *status);

int CallFunction(void *function,           ///< the function to call
                 cFunctionData *attribute, ///< attributes of the function
                 void *data) ///< private data for CCTK_CallFunction
{
    std::cout << "Calling: " << attribute->thorn << "::" << attribute->routine << " {" << std::endl;
    int const res = CCTK_CallFunction(function, attribute, data);
    std::cout << "}" << std::endl;
    return res;
}

extern "C"
void evol_here(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS_evol_here;
    std::cout << "evol here "<< std::endl;
}

void ScheduleTraverse(char const *const where, char const *const name,
                      cGH *const cctkGH) {
  CCTK_ScheduleTraverse(name, cctkGH, CallFunction);
}

int Evolve(tFleshConfig *const fc) 
{
    int conv_level = 0;
    cGH *cctkGH = fc->GH[conv_level];
    std::cout << "Evolve" << std::endl;
    ScheduleTraverse(where, "CCTK_PREREGRID", cctkGH);
    // do regrid
    ScheduleTraverse(where, "CCTK_POSTREGRID", cctkGH);
    // rotate time levels
    // t -> t + dt
    ScheduleTraverse(where, "CCTK_PRESTEP", cctkGH);
    ScheduleTraverse(where, "CCTK_EVOL", cctkGH);
    // Evolve finer grids recursively
    // Restrict from finer grids
    ScheduleTraverse(where, "CCTK_POSTRESTRICT", cctkGH);
    ScheduleTraverse(where, "CCTK_POSTSTEP", cctkGH);
    ScheduleTraverse(where, "CCTK_CHECKPOINT", cctkGH);
    // OutputGH
    ScheduleTraverse(where, "CCTK_ANALYSIS", cctkGH);
    return 0;
}

cGH *InstantiateGH(tFleshConfig *const fc,int conv_level)
{
    DECLARE_CCTK_PARAMETERS;

    cGH *thisGH = new cGH();
    std::cout << "cGH=" << thisGH << std::endl;
    thisGH->cctk_dim = CCTK_MaxGFDim();
    CCTKi_SetupGHExtensions(fc, conv_level, thisGH);

    /* Need this to be at least one otherwise the memory allocation will fail */
    int cctk_dim = thisGH->cctk_dim;
    if (thisGH->cctk_dim == 0)
    {
      cctk_dim = 1;
    }
    thisGH->cctk_iteration    = 0;
    thisGH->cctk_gsh          = new int[cctk_dim];
    thisGH->cctk_lsh          = new int[cctk_dim];
    thisGH->cctk_lbnd         = new int[cctk_dim];
    thisGH->cctk_ubnd         = new int[cctk_dim];
    thisGH->cctk_tile_min     = new int[cctk_dim];
    thisGH->cctk_tile_max     = new int[cctk_dim];

    thisGH->cctk_ash          = new int[cctk_dim];
    thisGH->cctk_to           = new int[cctk_dim];
    thisGH->cctk_from         = new int[cctk_dim];
    thisGH->cctk_bbox         = new int[2*cctk_dim];
    thisGH->cctk_nghostzones  = new int[2*cctk_dim];
    thisGH->cctk_levfac       = new int[cctk_dim];
    thisGH->cctk_levoff       = new int[cctk_dim];
    thisGH->cctk_levoffdenom  = new int[cctk_dim];
    thisGH->cctk_delta_space  = new CCTK_REAL8[cctk_dim];
    /* FIXME : Next line goes when coords are done properly */
    thisGH->cctk_origin_space = new CCTK_REAL8[cctk_dim];

    thisGH->cctk_delta_time = 1;
    thisGH->cctk_timefac = 1;
    thisGH->cctk_convlevel = 0;
    thisGH->data = new void**[CCTK_NumVars()];

    assert(not CCTK_IsFunctionAliased("MultiPatch_GetDomainSpecification"));
    //assert(CCTK_IsFunctionAliased("CoordBase_GetDomainSpecification"));

    // Ensure that CartGrid3D::type = "coordbase"
    if (CCTK_IsThornActive("CartGrid3D")) {
        int type;
        void const *const ptr = CCTK_ParameterGet("type", "CartGrid3D", &type);
        assert(ptr != 0);
        assert(type == PARAMETER_KEYWORD);
        char const *const coordtype = *static_cast<char const *const *>(ptr);
        if (not CCTK_EQUALS(coordtype, "coordbase")) {
            CCTK_ERROR("CartGrid3D is active. Please also set "
                       "CartGrid3D::type = \"coordbase\"");
        }
    }

    typedef double real3d[3];

    real3d physical_min, physical_max;
    real3d interior_min, interior_max;
    real3d exterior_min, exterior_max;
    real3d base_spacing;

    int ierr = GetDomainSpecification(
        cctk_dim, &physical_min[0], &physical_max[0], &interior_min[0],
        &interior_max[0], &exterior_min[0], &exterior_max[0], &base_spacing[0]);
    assert(not ierr);

    for(int i=0;i < cctk_dim;i++) {
        thisGH->cctk_ash[i] =
            std::ceil((exterior_max[i] - exterior_min[i])/base_spacing[i]);
        thisGH->cctk_lsh[i] =
            std::ceil((interior_max[i] - interior_min[i])/base_spacing[i]);
    }

    return thisGH;
}

void *MySetupGH(tFleshConfig *fc, int convLevel, cGH *cgh)
{
    std::cout << "MySetupGH was called" << std::endl;
    MyDriverData *md = new MyDriverData();
    std::cout << "md=" << md << std::endl;
    return md;
}

int MyInitGH(cGH *cgh)
{
    std::cout << "MyInitGH was called" << std::endl;
    return 0;
}

int Initialise(tFleshConfig *const fc) 
{
    DECLARE_CCTK_PARAMETERS;
    std::cout << "Initialise" << std::endl;

    int conv_level = 0;

    // We must initialize the cGH here. Borrow the default routine.
    cGH *cctkGH = InstantiateGH(fc, conv_level);
    fc->GH = new cGH*[1];
    fc->GH[0] = cctkGH;

    /* Initialise time */
    cctkGH->cctk_time = *(const CCTK_REAL *)
                  CCTK_ParameterGet ("cctk_initial_time", "Cactus", NULL);

    /* Initialise iteration number */
    cctkGH->cctk_iteration = 0;

    CCTKi_ScheduleGHInit(cctkGH);
    //CCTK_RegisterGHExtensionInitGH(handle, ghFunc);
    CCTKi_InitGHExtensions(cctkGH);
    MyDriverData *md = (MyDriverData *)CCTK_GHExtension(cctkGH, ExtensionName);

    ScheduleTraverse(where, "CCTK_WRAGH", cctkGH);
    ScheduleTraverse(where, "CCTK_PARAMCHECK", cctkGH);
    if(fc->recovered) {
        ScheduleTraverse(where, "CCTK_BASEGRID", cctkGH);
        ScheduleTraverse(where, "CCTK_RECOVER_VARIABLES", cctkGH);
        ScheduleTraverse(where, "CCTK_POST_RECOVER_VARIABLES", cctkGH);
    } else {
        ScheduleTraverse(where, "CCTK_PREREGRIDINITIAL", cctkGH);
        // setup grid hierarchy
        // Post
        ScheduleTraverse(where, "CCTK_PREREGRIDINITIAL", cctkGH);
        ScheduleTraverse(where, "CCTK_BASEGRID", cctkGH);
        ScheduleTraverse(where, "CCTK_INITIAL", cctkGH);
        ScheduleTraverse(where, "CCTK_POSTINITIAL", cctkGH);
        // initialize finer grids recursively
        ScheduleTraverse(where, "CCTK_POSTRESTRICTINITIAL", cctkGH);
        ScheduleTraverse(where, "CCTK_POSTPOSTINITIAL", cctkGH);
        ScheduleTraverse(where, "CCTK_POSTSTEP", cctkGH);
        ScheduleTraverse(where, "CCTK_CPINITIAL", cctkGH);
    }
    ScheduleTraverse(where, "CCTK_ANALYSIS", cctkGH);
    // OutputGH happens here
    return 0;
}

int Shutdown(tFleshConfig *const fc) 
{
    std::cout << "Shutdown" << std::endl;
    int conv_level = 0;
    cGH *cctkGH = fc->GH[conv_level];
    ScheduleTraverse(where, "CCTK_TERMINATE", cctkGH);
    //free(&cctkGH[0]); <-- wrong, causes double free
    delete cctkGH;
    ScheduleTraverse(where, "CCTK_SHUTDOWN", nullptr);
    return 0;
}
int GroupStorageDecrease(const cGH *cctkGH, int n_groups, const int *groups,
                       const int *requested_tls, int *status) {
    std::cout << "Calling Decrease...\n";
    return 0;
}
int GroupStorageIncrease(const cGH *cctkGH, int n_groups, const int *groups,
                       const int *requested_tls, int *status) {
    assert(cctkGH->data != nullptr);
    std::cout << "Calling De/Increase...\n";
    for(int i=0;i<n_groups;i++) {
        int ntls = requested_tls[i];
        int const group = groups[i];
        std::cout << "Increase for group " << group << std::endl;
        int const declared_tls = CCTK_DeclaredTimeLevelsGI(group);
        std::cout << "declared tls: " << declared_tls << std::endl;
        if(ntls < declared_tls)
            ntls = declared_tls;
        std::cout << "ntls: " << ntls << std::endl;

        int otls = 0;
        bool init = false;
        MyDriverData *md = (MyDriverData *)CCTK_GHExtension(cctkGH, ExtensionName);
        std::cout << "md: " << ((void*)md) << std::endl;
        assert(md != nullptr);
        auto fgroup = md->time_levels.find(group);
        if(fgroup == md->time_levels.end()) {
            init = true;
            std::cout << "init" << std::endl;
        } else {
            otls = fgroup->second;
            std::cout << "otls: " << otls << std::endl;
        }
        if(otls == ntls) {
            std::cout << "Nothing to do!" << std::endl;
            return 0;
        }
        md->time_levels[group] = ntls;

        cGroup gp;
        int ierr = CCTK_GroupData(group, &gp);
        int const firstvarindex = CCTK_FirstVarIndexI(group);
        for(int v=0;v<gp.numvars;v++) {
            int varindex = v + firstvarindex;
            std::cout << "varindex=" << varindex << std::endl;
            for(int j=0;j<ntls;j++) {
                std::cout << " allocate " << ntls << std::endl;
                void **new_data = new void*[ntls];
                int m = std::min(otls,ntls);
                for(int i=0;i<m;i++) {
                    std::cout << "  copy " << i << std::endl;
                    new_data[i] = cctkGH->data[varindex][i];
                }
                for(int i=m;i<ntls;i++) {
                    int dims = cctkGH->cctk_dim;
                    unsigned int sz = cctkGH->cctk_ash[dims-1];
                    sz = 1;
                    int gtype = CCTK_GroupTypeI(group);
                    int vtype = CCTK_VarTypeI(varindex);
                    int gsize = CCTK_VarTypeSize(vtype);
                    if(gsize <= 0) {
                        std::cout << "gtype=" << gtype << std::endl;
                        std::cout << "gsize=" << gsize << std::endl;
                        std::cout << "name=" << CCTK_FullGroupName(group) << std::endl;
                    }
                    assert(gsize > 0);
                    sz *= gsize;
                    if(gtype == CCTK_GF || gtype == CCTK_ARRAY) {
                        for(int dim=0;dim < dims;dim++) {
                            assert(cctkGH->cctk_ash[dim] > 0);
                            sz *= cctkGH->cctk_ash[dim];
                        }
                    }
                    while((sz % sizeof(long)) != 0)
                        sz++;
                    new_data[i] = new long[sz/sizeof(long)];
                }
                for(int i=m;i<otls;i++)
                    delete[]  ((long*)cctkGH->data[varindex][i]);
                if(!init) {
                    std::cout << "  delete" << std::endl;
                    delete[] cctkGH->data[varindex];
                }
                cctkGH->data[varindex] = new_data;
            }
        }
    }
    return 0;
}

int QueryGroupStorageB(const cGH *cctkGH, int group, const char *groupname) {
    std::cout << "Query storage" << std::endl;
    return 0;
}

extern "C"
int Driver_Startup()
{
  CCTK_RegisterBanner("This is a Custom Driver");

  CCTK_OverloadEvolve(Evolve);
  CCTK_OverloadInitialise(Initialise);
  CCTK_OverloadShutdown(Shutdown);

  CCTK_OverloadGroupStorageIncrease(GroupStorageIncrease);
  CCTK_OverloadGroupStorageDecrease(GroupStorageIncrease);

  CCTK_OverloadQueryGroupStorageB(QueryGroupStorageB);

  int GHExtension = CCTK_RegisterGHExtension(ExtensionName);
  CCTK_RegisterGHExtensionSetupGH(GHExtension, MySetupGH);
  CCTK_RegisterGHExtensionInitGH(GHExtension, MyInitGH);
  std::cout << "Registered Extension: " << ExtensionName << std::endl;

  #if 0
  CCTK_OverloadOutputGH(OutputGH);

  CCTK_OverloadSyncGroupsByDirI(SyncGroupsByDirI);

  // Don't need these two, older, call the next two
  CCTK_OverloadEnableGroupStorage(EnableGroupStorage);
  CCTK_OverloadDisableGroupStorage(DisableGroupStorage);

  // Need these two
  CCTK_OverloadGroupStorageIncrease(GroupStorageIncrease);
  CCTK_OverloadGroupStorageDecrease(GroupStorageDecrease);

  // Max number that were ever active
  CCTK_OverloadQueryMaxTimeLevels(QueryMaxTimeLevels);

  // Not really used
  CCTK_OverloadEnableGroupComm(EnableGroupComm);
  CCTK_OverloadDisableGroupComm(DisableGroupComm);

  CCTK_OverloadBarrier(Barrier);
  // CCTK_OverloadNamedBarrier (NamedBarrier);

  CCTK_OverloadExit((int (*)(cGH *, int))Exit);
  CCTK_OverloadAbort((int (*)(cGH *, int))Abort);

  CCTK_OverloadMyProc(MyProc);
  CCTK_OverloadnProcs(nProcs);

  // Return the size of a grid array in a particular format for Fortran
  CCTK_OverloadArrayGroupSizeB(ArrayGroupSizeB);
  // Whether a group has storage or not?
  CCTK_OverloadQueryGroupStorageB(QueryGroupStorageB);

  // Fill in lsh, etc. when grids change size
  // Sam and Roland did this for CarpetX
  CCTK_OverloadGroupDynamicData(GroupDynamicData);

  // Global constant for the number of dimensions? Global parameter?
  #endif

  return 0;
}
