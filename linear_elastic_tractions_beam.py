"""
Linear elasticity with pressure traction load on a surface and
constrained to one-dimensional motion.
"""
import numpy as nm


def linear_tension(ts, coor, mode=None, region=None, ig=None):
    if mode == 'qp':
    	print coor.shape[0] ;
    	print coor;
        #val = nm.tile( -1.0, (coor.shape[0], 1, 1))
        val = nm.zeros([coor.shape[0],3,1])
        #val[:] = [-1.0,0.,0.]
    	
    	val[:,1,:]=1.
    	print val;
    	print val.shape ;
        return {'val' : val}


# Define the function post_process, that will be called after the problem is
# solved.
def post_process(out, problem, state, extend=False):
    """
    This will be called after the problem is solved.

    Parameters
    ----------
    out : dict
        The output dictionary, where this function will store additional data.
    problem : ProblemDefinition instance
        The current ProblemDefinition instance.
    state : array
        The computed state vector, containing FE coefficients of all the
        unknown variables.
    extend : bool
        The flag indicating whether to extend the output data to the whole
        domain. It can be ignored if the problem is solved on the whole domain
        already.

    Returns
    -------
    out : dict
        The updated output dictionary. 
    """
    from sfepy.base.base import Struct

    # Cauchy strain averaged in elements.
    strain = problem.evaluate('de_cauchy_strain.i1.Omega( u )')
    out['cauchy_strain'] = Struct(name='output_data',
                                  mode='cell', data=strain,
                                  dofs=None)
    # Cauchy stress averaged in elements.
    stress = problem.evaluate('de_cauchy_stress.i1.Omega( solid.D, u )')
    out['cauchy_stress'] = Struct(name='output_data',
                                  mode='cell', data=stress,
                                  dofs=None)
    
    return out
     
    
    

def define():
    """Define the problem to solve."""
    from sfepy import data_dir

    filename_mesh = data_dir + '/meshes/3d/block.mesh'

    options = {
    	'post_process_hook' : 'post_process',
        'nls' : 'newton',
        'ls' : 'ls',
    }

    functions = {
        'linear_tension' : (linear_tension,),
    }

    fields = {
        'displacement': ('real', 3, 'Omega', 1),
    }

    materials = {
        'solid' : ({
            'lam' : 5.769,
            'mu' : 3.846,
        },),
        'load' : (None, 'linear_tension')
    }
	
    from sfepy.mechanics.matcoefs import stiffness_tensor_lame
	
    solid = materials['solid'][0]
    lam, mu = solid['lam'], solid['mu']
    solid.update({
		'D' : stiffness_tensor_lame(3, lam=lam, mu=mu),
	})    

    variables = {
        'u' : ('unknown field', 'displacement', 0),
        'v' : ('test field', 'displacement', 'u'),
    }

    regions = {
        'Omega' : ('all', {}),
        'Left' : ('nodes in (x < -4.99)', {}),
        'Right' : ('nodes in (x > 4.99)', {}),
    }

    ebcs = {
        'fixb' : ('Left', {'u.all' : 0.0}),
#        'fixt' : ('Right', {'u.[1,2]' : 0.0}),
    }
	#! Integrals
	#! ---------
	#! Define the integral type Volume/Surface and quadrature rule
	#! (here: dim=3, order=1).
    integrals = {
	    'i1' : ('v', 'gauss_o1_d3'),
	}
    ##
    # Balance of forces.
    equations = {
        'elasticity' :
        """dw_lin_elastic_iso.i1.Omega( solid.lam, solid.mu, v, u )
         = - dw_surface_ltr.i1.Right( load.val, v )""",
    }

    ##
    # Solvers etc.
    solvers = {
        'ls' : ('ls.scipy_direct', {}),
        'newton' : ('nls.newton',
                    { 'i_max'      : 1,
                      'eps_a'      : 1e-10,
                      'eps_r'      : 1.0,
                      'macheps'   : 1e-16,
                      # Linear system error < (eps_a * lin_red).
                      'lin_red'    : 1e-2,                
                      'ls_red'     : 0.1,
                      'ls_red_warp' : 0.001,
                      'ls_on'      : 1.1,
                      'ls_min'     : 1e-5,
                      'check'     : 0,
                      'delta'     : 1e-6,
                      'is_plot'    : False,
                      # 'nonlinear' or 'linear' (ignore i_max)
                      'problem'   : 'nonlinear'}),
    }

    ##
    # FE assembling parameters.
    fe = {
        'chunk_size' : 1000,
        'cache_override' : False,
    }

    return locals()
