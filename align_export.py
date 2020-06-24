import ast
from pymol import cmd, selector

"""
-fetch the structures(arg) 
-parseout the chains of interest(arg) 
-align in the same coordinate frame
-save as a pdb
"""
def chain_align_save(*args, **kwargs):
    args = [ast.literal_eval(kvpair) for kvpair in args]
    for pair in args:
        cmd.fetch(str.lower(pair[0]))
        create_subchain_object(pair[0], pair[1])
        cmd.delete(str.lower(pair[0]))
    
    for mobile in [ '{}.{}'.format(model, chain) for (model, chain) in args ][1:]:
        # cmd.align(mobile, "{}.{}".format(args[0][0], args[0][1]))
        cmd.super(mobile, "{}.{}".format(args[0][0], args[0][1]),reset=1,transform=1,quiet=0)

    cmd.reset()

    cmd.save('alignment.cif')
    

    

def create_subchain_object(pdbid, subchain):
    tempname = 'selection_{}.{}'.format(pdbid, subchain)
    cmd.select(tempname ,'m. {} and c. {}'.format( pdbid, subchain))
    cmd.create('{}.{}'.format(pdbid,subchain), tempname)
    cmd.delete(tempname)


cmd.extend('create_subchain_object', create_subchain_object)
cmd.extend('chain_align_save', chain_align_save);


