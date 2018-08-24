

def get_version_str(major,minor,minor_minor):
    return "%d.%d.%d"%(major,minor,minor_minor)

def get_output_loc(param):
    version_major   = param.get_int("version_major")
    version_minor   = param.get_int('version_minor')
    version_minor_minor = param.get_int('version_minor_minor')
    version_str = get_version_str(version_major,version_minor,version_minor_minor)
    output = param.get_string("output")
    return output.replace('${version}',version_str)
