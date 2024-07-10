import os
if os.name == 'nt':
	from .Release.pyrtklib5 import *
else:
	from .pyrtklib5 import *