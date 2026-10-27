python_blas <- function(num_threads) {
  reticulate::py_run_string(paste0("
from os import environ
environ['OMP_NUM_THREADS'] = '", as.integer(num_threads), "'
environ['OPENBLAS_NUM_THREADS'] = '", as.integer(num_threads), "'
"))
}


shutdown_joblib <- function(kill_workers=FALSE) {
  
  if(!reticulate::py_available(initialize=FALSE)){
    return(invisible(NULL))
  }
  
  kill_py <- if(isTRUE(kill_workers)) "True" else "False"
  
  reticulate::py_run_string(sprintf("
from joblib.externals.loky import reusable_executor

storage = getattr(reusable_executor, '_executor_storage', None)

executor = None
if storage is not None:
    executor = getattr(storage, 'executor', None)

if executor is not None:
    if hasattr(executor, 'terminate'):
        executor.terminate(kill_workers=%s)
    else:
        executor.shutdown(wait=True)

    # Drop the cached Python reference as well.
    storage.executor = None

    if hasattr(storage, 'executor_kwargs'):
        storage.executor_kwargs = None
", kill_py))
  
  invisible(NULL)
}
  

tn_python <- function(inner_threads=1L) {
  reticulate::py_run_string(paste0("
from joblib import parallel_config
import truncated_normal.truncated_normal as _tn

def tn_test(*args, **kwargs):
  # Limit BLAS/OpenMP/etc. in every loky subprocess.
  with parallel_config(
    backend='loky',
    inner_max_num_threads=", inner_threads, ",
    idle_worker_timeout=30
  ):
    return _tn.tn_test(*args, **kwargs)
"))
}


