"""
A set of common convinience functions
"""
#from IPython.testing.tools import full_path
import numpy as np
import os

# -----------------------------------------------------------------------------
def common_elements(list1, list2):
# -----------------------------------------------------------------------------
    """
    Calculates the number of common elements of 2 lists.
    """
    result = 0;
    for element in list1:
        if list2.__contains__(element):
            result += 1
    return result

# -----------------------------------------------------------------------------
def random_choice(options, weights = None, np_rnd = None):
# -----------------------------------------------------------------------------
    """
    Returns a random selection from the list od options. The probabilities of
    each option can be weighted. len(weights) must be the same as len(options).

    Not the fastest way, but it does the trick.
    """
    from bisect import bisect_right
    if weights:
        assert len(options) == len(weights)
    else:
        weights = [1]*len(options)

    if len(options) == 1:
        return options[0]

    if not np_rnd:
        np_rnd = np.random

    borders = [0]
    for weight in weights:
        borders.append(weight+borders[-1])
    borders = 1.*np.array(borders)/borders[-1]

    rnd_nr = np_rnd.random()

    # get to the correct bin with bisection
    return options[bisect_right(borders,rnd_nr)-1]

# -----------------------------------------------------------------------------
def get_default_function_parameters(func):
# -----------------------------------------------------------------------------
    '''
    Returns the default parameters of a given functions in a dictionary form.
    '''
    import inspect
    args, _, _, defaults = inspect.getargspec(func)
    n = len(defaults)
    result = {}
    for i,default in enumerate(defaults):
        result[args[len(args)-n+i]]=default
    return result

# -----------------------------------------------------------------------------
def save_figure_safe(full_path):
# -----------------------------------------------------------------------------
    """
    Makes sure the saving of the figure using pylab doesn't overwrite existing
    """
    import matplotlib.pylab as pyl

    if full_path.__contains__("."):
        file_type = full_path.split(".")[-1]
        name = full_path.rstrip("."+file_type)
    else:
        assert False, "ERROR: Specify the file type ..."

    path_candidate = full_path
    suffix=0
    while os.path.exists(path_candidate):
        path_candidate = "{0}_{1}.{2}".format(name, suffix, file_type)
        suffix += 1

    if path_candidate != full_path:
        print "WARNING: the file {0} already exists.".format(full_path)
        print "Saving the figure in the location {0}.".format(path_candidate)
    pyl.savefig(path_candidate)

# -----------------------------------------------------------------------------
def save_np_results_safe(full_path, np_array):
# -----------------------------------------------------------------------------
    """
    Makes sure the saving of the numpy results using pylab doesn't overwrite existing
    """
    if full_path.__contains__("."):
        file_type = full_path.split(".")[-1]
        name = full_path.rstrip("."+file_type)
    else:
        name = full_path
        file_type = "npy"

    path_candidate = full_path
    suffix=0
    while os.path.exists(path_candidate):
        path_candidate = "{0}_{1}.{2}".format(name, suffix, file_type)
        suffix += 1

    if path_candidate != full_path:
        print "WARNING: the file {0} already exists.".format(full_path)
        print "Saving the results in the location {0}.".format(path_candidate)
    np.save(path_candidate, np_array)

# -----------------------------------------------------------------------------
def save_pickle_safe(full_path, result):
# -----------------------------------------------------------------------------
    """
    Makes sure the pickling of the results using pylab doesn't overwrite
    existing files.
    Returns the full path used for saving.
    """
    import cPickle as pickle

    file_type = full_path.split(".")[-1]
    name = full_path.rstrip("."+file_type)

    path_candidate = full_path
    suffix=0
    while os.path.exists(path_candidate):
        path_candidate = "{0}_{1}.{2}".format(name, suffix, file_type)
        suffix += 1

    if path_candidate != full_path:
        print "WARNING: the file {0} already exists.".format(full_path)
        print "Pickling the results in the location {0}.".format(path_candidate)
    output=open(path_candidate, 'wb')
    pickle.dump(result, output,2)
    output.close()
    return path_candidate

# -----------------------------------------------------------------------------
def save_pickle(full_path, result):
# -----------------------------------------------------------------------------
    """
    Saving a pickle overwriting the exisitng files. Convinience method really.
    """
    import cPickle as pickle

    output=open(full_path, 'wb')
    pickle.dump(result, output,2)
    output.close()

# -----------------------------------------------------------------------------
def save_pickle_tmp(result):
# -----------------------------------------------------------------------------
    """
    Temporarly saving a pickle and returning the resulting path.
    """
    PATH = "/tmp/neurovivo_tmp.pkl"
    full_path = save_pickle_safe(PATH, result)
    return full_path

# -----------------------------------------------------------------------------
def mkdir_safe(path):
# -----------------------------------------------------------------------------
    '''
    Creates a folder. If the folder already exists it creates a new unique one.
    It returns the full path of the created folder.
    Does not work remotely.
    TODO: implement the remote version as well.
    '''
    if ('@' in path and ":" in path):
        assert False, "Not implemented for remote."
    prefix, folder_name = path.rsplit("/",1)
    path_candidate = path
    suffix = 0
    while os.path.exists(path_candidate):
        path_candidate = prefix + "/" + folder_name + str(suffix)
        suffix += 1
    mkdir_p(path_candidate)
    return path_candidate

# -----------------------------------------------------------------------------
def mkdir_p(path):
# -----------------------------------------------------------------------------
    '''
    Creates a folder or not if there is a folder already. sth like mkdir -p
    '''
    if not ('@' in path and ":" in path):
        try:
            os.makedirs(path)
        except OSError, e:
            if e.errno != 17:
                raise
    else:
        user_name=path.split("@")[0]
        file_name=path.split(":")[1]
        server_name=path.split("@")[1].split(":")[0]
        try:
            os.system("ssh {0}@{1} mkdir -p {2}".format(user_name, server_name, file_name))
        except OSError, e:
            if e.errno != 17:
                raise

# -----------------------------------------------------------------------------
def rm(path):
# -----------------------------------------------------------------------------
    '''
    tries to remove the path
    '''
    os.remove(path)

# -----------------------------------------------------------------------------
def rm_rf(path):
# -----------------------------------------------------------------------------
    '''
    tries to remove the path (does rm -rf - like command) and does not complain
    if the path does not exist.
    '''
    import shutil
    try:
        shutil.rmtree(path)
    except OSError:
        pass

# -----------------------------------------------------------------------------
def list_files(PATH, hidden=False, folders=False):
# -----------------------------------------------------------------------------
    '''
    returns a list of file names at the particular path.
    If hidden is True, the files starting with . are included.
    If folders is True, the folders are included.

    If the path is remote hidden is False and folders is True regardless what is
    specified.
    TODO: fix this a bit - at least add a warning
    '''
    # PATH has to end with / - even if there are two / it's ok ...
    PATH=PATH+"/"
    remote = False
    if not ('@' in PATH and ":" in PATH):
        files = os.listdir(PATH)
    else:
        remote = True
        user_name=PATH.split("@")[0]
        folder_path=PATH.split(":")[1]
        server_name=PATH.split("@")[1].split(":")[0]
        os.system("ssh {}@{} ls {} -1 >> tmp.txt".format(user_name, server_name, folder_path))
        f = open('tmp.txt', 'rb')
        files=[]
        for line in f:
            files.append(line[:-1])
        rm("tmp.txt")

    cp_files = files[:]
    if not hidden and not remote:
        for fajl in files:
            if fajl[0]=='.':
                cp_files.remove(fajl)

    if not folders and not remote:
        for fajl in files:
            if not os.path.isfile(PATH+fajl):
                cp_files.remove(fajl)

    #print cp_files
    return cp_files

# -----------------------------------------------------------------------------
def save_metadata(PATH, metadata):
# -----------------------------------------------------------------------------
    '''
    Saves the metadata. Used in the projects to save the relevant parameter
    values.
    '''
    mkdir_p(PATH + "/metadata/")
    save_pickle(PATH + "/metadata/metadata.pkl", metadata)

# -----------------------------------------------------------------------------
def load_metadata(PATH):
# -----------------------------------------------------------------------------
    '''
    Returns the metadata dictionary. It assumes the metadata is located in the
    file PATH/metadata/metadata.pkl which is where the save_metadata saves them
    to.
    '''
    import cPickle as pickle

    try:
        md = pickle.load(open(PATH+"/metadata/metadata.pkl"))
    except:
        raise
        assert False, "metadata missing. Expeceted metadata.pkl at {}".format(PATH+"/metadata/")

    return md

# -----------------------------------------------------------------------------
def load_pickle(PATH):
# -----------------------------------------------------------------------------
    """ Loads the pickled data. """
    import cPickle as pickle
    import tempfile

    if not ('@' in PATH and ":" in PATH):
        return pickle.load(open(PATH))
    else:
        file_name=PATH.split(":")[1].rsplit("/")[-1]
        tmp_path=tempfile.mkdtemp()
        os.system("scp {} {}".format(PATH, tmp_path))
        result = pickle.load(open(tmp_path+'/'+file_name))
        rm_rf(tmp_path)
        return result

# -----------------------------------------------------------------------------
def parameter_sweep(function, parameters, sweep_parameters={}, chosen_parameters={}, *args, **kwargs):
# -----------------------------------------------------------------------------
    """
    A recursive function that does the parameter sweep over the sweep parameters.
    PATH_RESULT - the absolute path to the folder where the results are saved
    function - the function used for the parameter sweep. Function should take
               a dictionary of parameters as an argument.
    parameters - a dictionary of parameters supplied to the function.
    sweep_parameters - a dictionary of parameters, where the key is a parameter
                       name and the value is a list of values for sweeping.
    chosen_parameters - internal parameter, used during recursion. No need to
                        bother about it.
    """
    sweep_dimension = len(sweep_parameters.keys())
    params=dict(parameters)
    chosen_params=dict(chosen_parameters)
    sweep_params=dict(sweep_parameters)
    if sweep_dimension > 0:
        still_open_keys = sweep_params.keys()
        still_open_keys.sort()
        exploded_param = still_open_keys[0]
        exploded_values = sweep_params.pop(exploded_param)
        exploded_values.sort()
        for value in exploded_values:
            params[exploded_param]=value
            chosen_params[exploded_param]=value # I realise this is storing more then needed, but helps with debugning
            parameter_sweep(function, params, sweep_params, chosen_params, *args, **kwargs)
    else:
        if "TEST" in parameters.keys():
            print "Starting with:"
            pretty_dict_print(chosen_params)
        params["chosen_sweep_params"]=chosen_params
        function(params, *args, **kwargs)

# -----------------------------------------------------------------------------
def save_with_info_string(chosen_params, result, PATH, name="result"):
# -----------------------------------------------------------------------------
    """
    Given a path for results it saves them as a pickled dictionary in a file
    with a formatted informative name based on the chosen_params.
    It returns the full file-name and path of the saved result.
    """
    info_string = create_info_string(chosen_params)
    save_pickle(PATH+"{0}{1}.pickle".format(name, info_string), result)
    return PATH+"{0}{1}.pickle".format(name, info_string)

# -----------------------------------------------------------------------------
def create_info_string(chosen_params):
# -----------------------------------------------------------------------------
    """
    Given a set of chosen parameters it creates an informative string with the
    info.
    """
    info_string=""
    keys = chosen_params.keys()
    keys.sort()
    for key in keys:
        info_string = "_".join([info_string,"{0}-{1}".format(key,chosen_params[key])])
    return info_string

# -----------------------------------------------------------------------------
def get_chosen_file_name(PATH_RESULTS, params):
# -----------------------------------------------------------------------------
    """
    Used when doing parameter sweeps.
    Given a PATH with files that represent results for different parameter
    values it returns a target file name.

    The first name that has all the parameters correctly specified is returned.
    No checks are done for cases where more files are correct.

    Looks very messy:
    TODO: clean it up
    """

    all_file_names=list_files(PATH_RESULTS)
    for file_name in all_file_names:
        got_one = True
        # cut out the part after the last "." that is the file extension
        file_name_split=file_name.rsplit(".",1)[0]
        parts = file_name_split.split("~")
        file_params={}
        for part in parts[1:]:
            [key, value] = part.split("-")
            try:
                value = float(value)
            except:
                pass
            file_params[key]=value
        for i in params.keys():
            # strings get treated differently then numbers
            if isinstance(params[i], str):
                if not params[i]==file_params[i]:
                    got_one=False
                    break
            else:
                if not round(params[i]-file_params[i],7)==0:
                    got_one=False
                    break
        # this means we had no mistakes in the parameters, so we return
        if got_one:
            return file_name

    assert False, "Did not find a resulting file name for parameters: {} in the path {}".format(params, PATH_RESULTS)

# -----------------------------------------------------------------------------
def load_results(PATH_RESULTS, chosen_params):
# -----------------------------------------------------------------------------
    """
    Given a path with files that represent results for different parameter
    values it loads and returns the result file (a dictionary).
    """
    file_name = get_chosen_file_name(PATH_RESULTS, chosen_params)
    return load_pickle(PATH_RESULTS+"/"+file_name)

# -----------------------------------------------------------------------------
def load_parameters(PATH_RESULTS, chosen_params):
# -----------------------------------------------------------------------------
    """
    Given a path with files that represent results for different parameter
    values it loads and returns the result file (a dictionary).
    """
    file_name = get_chosen_file_name(PATH_RESULTS + "/parameters/", chosen_params)
    return load_pickle(PATH_RESULTS+"/parameters/"+file_name)

# -----------------------------------------------------------------------------
def save_yaml(full_path, dictionary):
# -----------------------------------------------------------------------------
    """
    Saving a yaml file overwriting existing files. Convinience.
    """
    import yaml
    output=open(full_path, 'w')
    yaml.dump(dictionary, output)
    output.close()

# -----------------------------------------------------------------------------
def load_yaml(full_path, dictionary):
# -----------------------------------------------------------------------------
    """
    Loading a yaml dictionary. Convinience.
    """
    import yaml
    input_file=open(full_path, 'rb')
    yaml.load(dictionary, input_file)
    input_file.close()

# -----------------------------------------------------------------------------
def flatten(listOfLists):
# -----------------------------------------------------------------------------
    '''Flatten one level of nesting'''
    from itertools import chain
    return list(chain.from_iterable(listOfLists))

# -----------------------------------------------------------------------------
class Bunch(object):
# -----------------------------------------------------------------------------
    '''
    A helper class. Making a dictionary adict into an object. Once the object
    is created, each dictonary element can be accessed by object.element .
    '''
    def __init__(self, adict):
        self.__dict__.update(adict)

# -----------------------------------------------------------------------------
def pretty_dict_print(dictionary):
# -----------------------------------------------------------------------------
    '''
    A line by line print of a dictionary.
    '''
    for key in dictionary.keys():
        print "{}: {}".format(key, dictionary[key])

# -----------------------------------------------------------------------------
def make_report(PATH_RESULTS, name="", threshold=True):
# -----------------------------------------------------------------------------
    """
    Given the project folder with results it creates a singe PDF with the
    metadata and all the .pdf figures for that project.

    If threshold is True, only the files that contain the word "report" will be
    included.
    """
    PATH_REPORT = PATH_RESULTS+"/report/"
    mkdir_p(PATH_REPORT)
    files_for_report=[]
    metadata=None
    try:
        metadata=load_metadata(PATH_RESULTS)
    except:
        pass
    if not metadata is None:
        metadata_str=str(metadata)
        metadata_str2=metadata_str.replace(',', '<br/>')
        from reportlab.pdfgen import canvas
        from reportlab.lib.styles import getSampleStyleSheet
        from reportlab.platypus import Paragraph
        c = canvas.Canvas(PATH_RESULTS +"/metadata/metadata.pdf")
        styleSheet = getSampleStyleSheet()
        style = styleSheet['BodyText']
        P=Paragraph(metadata_str2,style)
        aW = 480    # available width and height
        aH = 100
        P.wrap(aW, aH)    # find required space
        P.drawOn(c,10,aH)
        c.showPage()
        c.save()
    from pyPdf import PdfFileWriter, PdfFileReader
    output = PdfFileWriter()
    try:
        figure_pdfs = list_files(PATH_RESULTS+"/figures/")
    except:
        figure_pdfs = []
    # we only really take the pdfs
    for fajl in figure_pdfs:
        if ".pdf" in fajl:
            files_for_report.append(PATH_RESULTS+"/figures/"+fajl)
            if threshold and not "report" in fajl:
                    files_for_report.pop(-1)
    files_for_report.sort()
    if not metadata is None:
        files_for_report_tmp=[PATH_RESULTS +"/metadata/metadata.pdf"]
    files_for_report_tmp.extend(files_for_report)

    for fajl in files_for_report_tmp:
        inputt = PdfFileReader(file(fajl, "rb"))
        for i in xrange(inputt.getNumPages()):
            output.addPage(inputt.getPage(i))
    outputStream = file(PATH_REPORT + "/report-{}.pdf".format(name), "wb")
    output.write(outputStream)
    outputStream.close()

# -----------------------------------------------------------------------------
def nextpow2(n):
# -----------------------------------------------------------------------------
    """
    Given a number n Returns the next bigger integer which is the power of 2.
    """
    m_f = np.log2(n)
    m_i = int(np.ceil(m_f))
    return 2**m_i


# -----------------------------------------------------------------------------
def send_email(email_text, subject="Simulation Finished"):
# -----------------------------------------------------------------------------
    """
    Sends me email. Used for notifications that the simulations have finished.
    """
    import smtplib
    from email.MIMEText import MIMEText
    from email.MIMEMultipart import MIMEMultipart

    gmailUser = 'mpelko.smtp@gmail.com'
    gmailPassword = 'smtpsmtp'
    recipient = 'mpelko@gmail.com'
    msg = MIMEMultipart()
    msg['From'] = gmailUser
    msg['To'] = recipient
    msg['Subject'] = subject
    msg.attach(MIMEText(email_text))
    mailServer = smtplib.SMTP('smtp.gmail.com', 587)
    mailServer.ehlo()
    mailServer.starttls()
    mailServer.ehlo()
    mailServer.login(gmailUser, gmailPassword)
    mailServer.sendmail(gmailUser, recipient, msg.as_string())
    mailServer.close()
    print('Sent email to %s' % recipient)
