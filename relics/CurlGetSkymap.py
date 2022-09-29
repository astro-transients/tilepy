# Download a LIGO-Virgo skymap using curl.
# Roy Williams <roy.williams@ligo.org>
 
# Add a line like this to your ~/.netrc file:
#    machine gracedb.ligo.org login albert.einstein@LIGO.ORG password JKJKJKJkjkjkjkj1223
# You can get the password from:
#    https://gracedb.ligo.org/gracedb/options/manage_password
 
import os
 
def download_skymap(grace_id, name='bayestar', base_url='https://gracedb.ligo.org/gracedb'):
    """
    Download a sky map from GraceDb.
 
    Parameters
    ----------
    grace_id : string
        The GraceDb identifier of the event
    name : string
        The filename, less the '.fits.gz' extension
    base_url : string
        The base URL of the GraceDb server. By default, the production server.
    """
 
    command = 'curl --netrc {base_url}/apibasic/events/{grace_id}/files/{name}.fits.gz -o {grace_id}_{name}.fits.gz'.format(
        base_url=base_url, grace_id=grace_id, name=name)
    print(command)
    os.system(command)
 
 
# Example using the function above:
if __name__ == '__main__':
    #grace_id = 'T125738'        # identifier for the event
    name='bayestar'
    grace_id = 'G268556'        # identifier for the event
    #name = 'lalinference'
    download_skymap(grace_id, name)   # fetch the sky map


#https://gracedb.ligo.org/events/G268556