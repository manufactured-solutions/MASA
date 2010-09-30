#define libmesh_write_traceout()   do { std::stringstream outname; outname << "traceout_" << libMesh::processor_id() << '_' << getpid() << ".txt"; std::ofstream traceout(outname.str().c_str()); libMesh::print_trace(traceout); } while(0)

#define libmesh_error()    do { libmesh_write_traceout(); libmesh_here(); LIBMESH_THROW(libMesh::LogicError()); } while(0)

#ifdef LIBMESH_ENABLE_EXCEPTIONS
#define LIBMESH_THROW(e) do { throw e; } while (0)
#else
#define LIBMESH_THROW(e) do { std::abort(); } while (0)

