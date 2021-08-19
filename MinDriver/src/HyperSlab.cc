#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <assert.h>

typedef CCTK_INT (*conversion_fn_ptr)(
    CCTK_INT const nelems, CCTK_INT const src_stride, CCTK_INT const dst_stride,
    CCTK_INT const src_type, CCTK_INT const dst_type,
    CCTK_POINTER_TO_CONST const from, CCTK_POINTER const to);

CCTK_INT
MinDriverSlab_GetList(CCTK_POINTER_TO_CONST const cctkGH_,
                   CCTK_INT const mapping_handle, CCTK_INT const num_arrays,
                   CCTK_INT const *const procs, CCTK_INT const *const vindices,
                   CCTK_INT const *const timelevels,
                   CCTK_INT const *const hdatatypes,
                   CCTK_POINTER const *const hdata, CCTK_INT *const retvals) {
  cGH const *const cctkGH = (cGH const *)cctkGH_;
  assert("Not implemented"==0);
  return 0;
}

CCTK_INT
MinDriverSlab_FreeMapping(CCTK_INT const mapping_handle) {
  assert("Not implemented"==0);
  return 0;
}

CCTK_INT
MinDriverSlab_GlobalMappingByIndex(
    CCTK_POINTER_TO_CONST const cctkGH_, CCTK_INT const vindex,
    CCTK_INT const hdim, CCTK_INT const *const direction,
    CCTK_INT const *const origin, CCTK_INT const *const extent,
    CCTK_INT const *const downsample_, CCTK_INT const table_handle,
    conversion_fn_ptr const conversion_fn, CCTK_INT *const hsize) {
  cGH const *const cctkGH = (cGH const *)cctkGH_;
  assert("Not implemented"==0);
  return 0;
}
