using Logging
using LoggingExtras
using Test
using TreeKnit

@testset "Loggers skip messages below their levels (verbosity $v)" for v in -1:3
	tee = TeeLogger(TreeKnit.get_loggers(v, IOBuffer())...)
	@test Logging.min_enabled_level(tee) == min(LogLevel(-1), LogLevel(-v))
end
