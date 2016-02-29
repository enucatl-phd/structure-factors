"""Default options for the logger."""

config_dictionary = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'standard': {
            'format':
            '%(asctime)s %(levelname)-5.5s [%(name)s] %(message)s'
            },
        },
    'handlers': {
        'default': {
            'level': 'INFO',
            'class': 'logging.StreamHandler',
            'formatter': 'standard',
            },
        },
    'loggers': {
        '': {
            'handlers': ['default'],
            'level': 'INFO',
            'propagate': True
            },
        }
    }
