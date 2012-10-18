#include <sstream>

#include "pan.h"

pan::internal::streambuf_to_android_log pan::android_sb;
std::ostream pan::lout(0);

void pan::init_log() 
{
    lout.rdbuf(&android_sb);
}

namespace pan
{
namespace internal
{
streambuf_to_android_log::streambuf_to_android_log ( std::size_t buff_sz, std::size_t put_back ) :
    put_back_ ( std::max ( put_back, size_t ( 1 ) ) ),
    buffer_ ( std::max ( buff_sz, put_back_ ) + put_back_ )
{
}
void streambuf_to_android_log::append ( int of, char* first, char* last )
{
    size_t size = std::distance ( first, last );

    //tmp_.clear();

    while ( first != last ) {

        if ( *first == '\n' ) {
            tmp_.push_back ( '\0' );
            __android_log_print ( ANDROID_LOG_INFO, "pan::log", tmp_.data() );
            tmp_.clear();
        } else {
            tmp_.push_back ( *first );
        }
        ++first;
    }

    // basically do the same once more for the overflow. TODO: get rid of the duplicate code somehow
    if ( of != 0 ) {
        if ( of == '\n' ) {
            tmp_.push_back ( '\0' );
            __android_log_print ( ANDROID_LOG_INFO, "pan::log", tmp_.data() );
            tmp_.clear();
        } else {
            tmp_.push_back ( of );
        }
    }
}
int streambuf_to_android_log::overflow ( int c )
{
    append ( c, pbase(), epptr() );

    setp ( &buffer_.front(), ( &buffer_.back() ) + 1 );

    return 1;

}
int streambuf_to_android_log::sync()
{
    append ( 0, pbase(), pptr() );
    setp ( &buffer_.front(), ( &buffer_.back() ) + 1 );
    return 0;

}
} // namespace pan::internal


std::string sl_error_exception::err_str ( SLresult res, const char* file, int line ) throw()
{
    std::stringstream ss;
    ss << "OpenSL error: " << errcode_string ( res ) << " at " << file << ":" << line;

    return ss.str();
}
const char* sl_error_exception::errcode_string ( SLresult res )
{

    switch ( res ) {

    case SL_RESULT_SUCCESS:
        return "SL_RESULT_SUCCESS";
    case SL_RESULT_PRECONDITIONS_VIOLATED:
        return "SL_RESULT_PRECONDITIONS_VIOLATED";
    case SL_RESULT_PARAMETER_INVALID:
        return "SL_RESULT_PARAMETER_INVALID";
    case SL_RESULT_MEMORY_FAILURE:
        return "SL_RESULT_MEMORY_FAILURE";
    case SL_RESULT_RESOURCE_ERROR:
        return "SL_RESULT_RESOURCE_ERROR";
    case SL_RESULT_RESOURCE_LOST:
        return "SL_RESULT_RESOURCE_LOST";
    case SL_RESULT_IO_ERROR:
        return "SL_RESULT_IO_ERROR";
    case SL_RESULT_BUFFER_INSUFFICIENT:
        return "SL_RESULT_BUFFER_INSUFFICIENT";
    case SL_RESULT_CONTENT_CORRUPTED:
        return "SL_RESULT_CONTENT_CORRUPTED";
    case SL_RESULT_CONTENT_UNSUPPORTED:
        return "SL_RESULT_CONTENT_UNSUPPORTED";
    case SL_RESULT_CONTENT_NOT_FOUND:
        return "SL_RESULT_CONTENT_NOT_FOUND";
    case SL_RESULT_PERMISSION_DENIED:
        return "SL_RESULT_PERMISSION_DENIED";
    case SL_RESULT_FEATURE_UNSUPPORTED:
        return "SL_RESULT_FEATURE_UNSUPPORTED";
    case SL_RESULT_INTERNAL_ERROR:
        return "SL_RESULT_INTERNAL_ERROR";
    case SL_RESULT_UNKNOWN_ERROR:
        return "SL_RESULT_UNKNOWN_ERROR";
    case SL_RESULT_OPERATION_ABORTED:
        return "SL_RESULT_OPERATION_ABORTED";
    case SL_RESULT_CONTROL_LOST:
        return "SL_RESULT_CONTROL_LOST";
    default:
        return "(unknown)";
    }
}
async_audio_output::async_audio_output()
 : have_fill_func_(false),
   have_fill_float_func_(false)
{
    
    const size_t buf_size = 32;
    
    for( auto &buf: bufs ) {
        buf.resize(buf_size);
    }
    buf_it = bufs.begin();
    
    float_buf.resize(buf_size);
    
    
    
    auto res = slCreateEngine ( &obj_engine, 0, 0, 0, 0, 0 );
    res = ( *obj_engine )->Realize ( obj_engine, SL_BOOLEAN_FALSE );
    check_sl_error ( res );
    res = ( *obj_engine )->GetInterface ( obj_engine, SL_IID_ENGINE, &if_engine );
    check_sl_error ( res );

    pan::lout << "SL_IID_ANDROIDSIMPLEBUFFERQUEUE " << SL_IID_ANDROIDSIMPLEBUFFERQUEUE << "\n";
    const SLInterfaceID ids[] = {SL_IID_VOLUME};
    const SLboolean req[] = {SL_BOOLEAN_FALSE};
    res = ( *if_engine )->CreateOutputMix ( if_engine, &obj_output_mix, 1, ids, req );
    check_sl_error ( res );
    res = ( *obj_output_mix )->Realize ( obj_output_mix, SL_BOOLEAN_FALSE );
    check_sl_error ( res );

    SLDataLocator_AndroidSimpleBufferQueue loc_buf_queue = {SL_DATALOCATOR_ANDROIDSIMPLEBUFFERQUEUE, 2};

    SLuint32 channels = 2;
    SLuint32 sample_rate = 44100 * 1000;
    SLuint32 speakers = 0;

    SLDataFormat_PCM format_pcm = {SL_DATAFORMAT_PCM, channels,sample_rate, SL_PCMSAMPLEFORMAT_FIXED_16, SL_PCMSAMPLEFORMAT_FIXED_16, speakers, SL_BYTEORDER_LITTLEENDIAN};
    SLDataSource audio_src = {&loc_buf_queue, &format_pcm};
    SLDataLocator_OutputMix loc_outmix = {SL_DATALOCATOR_OUTPUTMIX, obj_output_mix };
    SLDataSink audio_sink = {&loc_outmix, 0};

    const SLInterfaceID ids1[] = {SL_IID_ANDROIDSIMPLEBUFFERQUEUE};
    const SLboolean req1[] = {SL_BOOLEAN_TRUE};


    res = ( *if_engine )->CreateAudioPlayer ( if_engine, &obj_player, &audio_src, &audio_sink, 1, ids1, req1 );
    check_sl_error ( res );

    res = ( *obj_player )->Realize ( obj_player, SL_BOOLEAN_FALSE );
    check_sl_error ( res );


    res = ( *obj_player )->GetInterface ( obj_player, SL_IID_PLAY, &if_play );
    check_sl_error ( res );



    res = ( *obj_player )->GetInterface ( obj_player, SL_IID_ANDROIDSIMPLEBUFFERQUEUE, &if_bqueue );
    check_sl_error ( res );
    res = ( *if_bqueue )->RegisterCallback ( if_bqueue, bq_callback, this );
    check_sl_error ( res );

}

void async_audio_output::bq_callback ( SLAndroidSimpleBufferQueueItf caller, void* data )
{
    async_audio_output *self = reinterpret_cast<async_audio_output *> ( data );

    assert ( self != nullptr );
    assert ( caller == self->if_bqueue );
//     pan::lout << "callback" << std::endl;

    
    
    
    if( self->buf_it == self->bufs.end() ) {
        self->buf_it = self->bufs.begin();
    }
    
    auto &buf = *(self->buf_it);
    ++self->buf_it;
    
    if( self->have_fill_func_ ) {
        //assert ( self->fill_buffer_func_ != nullptr );
        ( self->fill_buffer_func_ ) ( &(*buf.begin()), &(*buf.end()) ); // FIXME: do checked STL impls like dereferencing ::end()?
        
        
    } else if( self->have_fill_float_func_ ) {
        ( self->fill_buffer_float_func_ ) ( &(*self->float_buf.begin()), &(*self->float_buf.end()) ); // FIXME: do checked STL impls like dereferencing ::end()?
        
        assert( buf.size() == self->float_buf.size() );
        std::transform( self->float_buf.begin(), self->float_buf.end(), buf.begin(), []( float f ) {
            return int16_t(f * 255 * 127);
            //return int16_t(f * 8000.0);
        } );
//         lout << "fill float" << std::endl;
    }
    
    
    ( *caller )->Enqueue ( caller, buf.data(), buf.size() * sizeof ( int16_t ) );
}
void async_audio_output::start()
{
    auto res = ( *if_play )->SetPlayState ( if_play, SL_PLAYSTATE_PLAYING );
    check_sl_error ( res );

    bq_callback ( if_bqueue, this );



}
async_audio_output::~async_audio_output()
{
    throw std::runtime_error ( "~async_audio_output(): out of order" );
}

} // namespace pan
