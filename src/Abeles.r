#include "XOPStandardHeaders.r"

resource 'vers' (1) {						/* XOP version info */
	0x01, 0x10, final, 0x00, 0,				/* version bytes and country integer */
	"1.10",
	"1.10, � 1996-2004 WaveMetrics, Inc., all rights reserved."
};

resource 'vers' (2) {						/* Igor version info */
	0x08, 0x00, release, 0x00, 0,			/* version bytes and country integer */
	"8.00",
	"(for Igor Pro 8.00 or later)"
};

resource 'STR#' (1100) {					/* custom error messages */
	{
	"Abeles requires Igor Pro 8.0 or later.",
	"Wave does not exist",
	"Coefficient wave must be single or double precision floating point",
	"The coefficient wave has the wrong number of parameters. Wavelength = 4*w[0]+6+4*mullayers",
	"Waves not same length",
	"Requires double precision wave",
	}
};

resource 'STR#' (1101) {					// Misc strings for XOP.
	{
		"-1",								// This item is no longer supported by the Carbon XOP Toolkit.
		"No Menu Item",						// This item is no longer supported by the Carbon XOP Toolkit.
		"Abeles Help",					// Name of XOP's help file.
	}
};

resource 'XOPI' (1100) {
	XOP_VERSION,							// XOP protocol version.
	DEV_SYS_CODE,							// Development system information.
	XOP_FEATURE_FLAGS,						// Obsolete - set to zero.
	XOPI_RESERVED,							// Obsolete - set to zero.
	XOP_TOOLKIT_VERSION,					// XOP Toolkit version.
};

resource 'XOPF' (1100) {
	{
		"Abeles",					// Function name.
		F_EXP | F_THREADSAFE | F_EXTERNAL,				// Function category,
		NT_FP64,
		{						// Return value type.
			NT_FP64 | WAVE_TYPE,		// Double precision wave (coefficient wave).
			NT_FP64,					// Double precision x.
		},
		
		"Abeles_imag",					// Function name.
		F_EXP | F_THREADSAFE | F_EXTERNAL,				// Function category,
		NT_FP64,	
		{					// Return value type.
			NT_FP64 | WAVE_TYPE,		// Double precision wave (coefficient wave).
			NT_FP64,					// Double precision x.
		},
			
		"AbelesAll",					// Function name.
		F_EXP | F_THREADSAFE | F_EXTERNAL,				// Function category,
		NT_FP64,						// Return value type.
		{
			NT_FP64 + WAVE_TYPE,		// Double precision wave (coefficient wave).
			NT_FP64 + WAVE_TYPE,		// Double precision wave (y wave).
			NT_FP64 + WAVE_TYPE,		// Double precision wave (x wave).
		},
		
		"Abeles_imagAll",					// Function name.
		F_EXP | F_THREADSAFE | F_EXTERNAL,				// Function category,
		NT_FP64,
		{						// Return value type.
			NT_FP64 + WAVE_TYPE,		// Double precision wave (coefficient wave).
			NT_FP64 + WAVE_TYPE,		// Double precision wave (y wave).
			NT_FP64 + WAVE_TYPE,		// Double precision wave (x wave).
		},

	}
};
