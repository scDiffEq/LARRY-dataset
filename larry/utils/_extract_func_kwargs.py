import inspect


def _extract_func_params(func):
    return list(inspect.signature(func).parameters.keys())


def extract_func_kwargs(func, kwargs: dict, ignore: list=[]):
    
    """
    Parameters:
    -----------
    func
        type: Any
    
    kwargs
        type: dict
    
    ignore
        type: list
        default: []
        
    Returns:
    --------
    func_kwargs
        type: dict
    """
    
    func_kwargs = {}
    params = _extract_func_params(func)
    
    for k, v in kwargs.items():
        if (k in params) and (not k in ignore):
            func_kwargs[k] = v
            
    return func_kwargs