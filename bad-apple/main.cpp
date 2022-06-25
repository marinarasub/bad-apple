#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include <stdexcept>
#include <thread>
#include <chrono>
#include <cassert>
#include <thread>
#include <mutex>

#include "png.h"
#include "FLAC++/all.h"
#include "portaudio.h"

constexpr size_t PNG_HEADER_BYTES = 8;

typedef unsigned char byte;

uint32_t uint_max_val(uint8_t bit_depth)
{
    uint32_t max_int = 0;
    for (uint8_t i = 0; i < bit_depth; i++)
    {
        max_int = (max_int << 1) | 1;
    }
    return max_int;
}

struct audio_sample
{
    float left; /**< Left channel sample. */
    float right; /**< Right channel sample. */

    audio_sample operator+(const audio_sample& other) const
    {
        audio_sample s{};
        s.left = this->left + other.left;
        s.right = this->right + other.right;
        s.clamp();
        return s;
    }

    void clamp();
};


class AudioTrack
{
    static constexpr byte MONO = 1;
    static constexpr byte STEREO = 2;
    static constexpr unsigned int DEFAULT_SAMPLE_RATE = 44100;

public:

    AudioTrack();
    AudioTrack(unsigned long length);
    AudioTrack(unsigned long length, unsigned long sample_rate); // in samples
    AudioTrack(const AudioTrack& other);

    ~AudioTrack();

    AudioTrack operator=(const AudioTrack& other);
    AudioTrack operator+(const AudioTrack& other) const;

    unsigned int sampleRate() const;
    unsigned long sampleLength() const;
    byte channels() const;
    void resetPosition();
    unsigned long currentPos() const;
    bool ended() const;
    audio_sample nextSample();
    audio_sample getSample(unsigned long idx) const;
    void setSample(unsigned long idx, audio_sample data);
    void resample(unsigned long new_sample_rate);
    void merge(AudioTrack& track);

protected:
    void copy(const AudioTrack& other);
    void clear();
    audio_sample* audio; /**< Array of audio samples of the track */

private:

    void allocSamples(size_t length);

    unsigned int sample_rate; /**< Intended sample rate for playback*/
    unsigned long current_pos; /**< Cursor for current position in track */
    unsigned long length; /**< Length of audio sample array */
    byte num_tracks; /**< Number of audio channels in track */

};

constexpr int START_POS = 0;

// TODO variable sample rate
AudioTrack::AudioTrack()
{
    audio = nullptr;
    length = 0;
    num_tracks = 0;
    sample_rate = 0;
    current_pos = START_POS;
}

AudioTrack::AudioTrack(unsigned long length)
    : AudioTrack(length, DEFAULT_SAMPLE_RATE)
{
}

AudioTrack::AudioTrack(unsigned long length, unsigned long sample_rate)
    : num_tracks(STEREO), sample_rate(sample_rate), current_pos(START_POS)
{
    if (length == 0) throw std::invalid_argument("Length cannot be 0");
    allocSamples(length);
}

AudioTrack::AudioTrack(const AudioTrack& other) : AudioTrack()
{
    copy(other);
}

AudioTrack::~AudioTrack()
{
    clear();
}

AudioTrack AudioTrack::operator=(const AudioTrack& other)
{
    if (this != &other) copy(other);
    return *this;
}

AudioTrack AudioTrack::operator+(const AudioTrack& other) const
{
    unsigned long max_len = this->sampleLength() > other.sampleLength() ?
        this->sampleLength() : other.sampleLength();
    AudioTrack mix(max_len); // TODO pass in sample rate etc.
    for (unsigned long i = 0; i < max_len; i++)
    {
        audio_sample add = this->getSample(i) + other.getSample(i);
        mix.setSample(i, add);
    }
    return mix;
}

unsigned int AudioTrack::sampleRate() const
{
    return sample_rate;
}

unsigned long AudioTrack::sampleLength() const
{
    return length;
}

byte AudioTrack::channels() const
{
    return num_tracks;
}

void AudioTrack::resetPosition()
{
    current_pos = START_POS;
}

unsigned long AudioTrack::currentPos() const
{
    return current_pos;
}

bool AudioTrack::ended() const
{
    return current_pos >= length;
}

audio_sample AudioTrack::nextSample()
{
    audio_sample next{ 0, 0 };
    if (!ended()) {
        next = audio[currentPos()];
        current_pos++;
    }
    return next;
}

audio_sample AudioTrack::getSample(unsigned long idx) const
{
    if (idx >= sampleLength()) return audio_sample{ 0, 0 };
    return audio[idx];
}

void AudioTrack::setSample(unsigned long idx, audio_sample data)
{
    data.clamp();
    audio[idx] = data;
}

void AudioTrack::resample(unsigned long new_sample_rate)
{
    if (new_sample_rate == 0) throw std::invalid_argument("sample rate cannot be < 1");
    if (!audio) throw std::runtime_error("no audio data to resample");

    double ratio = (double)sample_rate / new_sample_rate;
    unsigned long new_len = (unsigned long)(1.0 / ratio * sampleLength());
    if (new_len == 0) throw std::runtime_error("error: resampled track length is 0");

    AudioTrack temp(new_len, new_sample_rate);
    for (unsigned long i = 0; i < new_len; i++)
    {
        // nearest neighbor resample
        temp.setSample(i, getSample((unsigned long)floor(ratio * i)));
    }
    *this = temp;
}

void AudioTrack::merge(AudioTrack& other)
{
    *this = *this + other;
}


void AudioTrack::copy(const AudioTrack& other)
{
    clear();
    allocSamples(other.length);
    this->num_tracks = other.num_tracks;
    this->sample_rate = other.sampleRate();
    for (unsigned long i = 0; i < length; i++)
    {
        this->setSample(i, other.getSample(i));
    }
}

void AudioTrack::clear()
{
    if (audio != nullptr)
    {
        delete[] audio;
        audio = nullptr;
    }
}

void AudioTrack::allocSamples(size_t length)
{
    this->length = length;
    audio = new audio_sample[length];
}

void audio_sample::clamp()
{
    left = fmin(1.0, fmax(-1.0, left));
    right = fmin(1.0, fmax(-1.0, right));
}


namespace FLACAPI = FLAC;

class FLACTrack : public AudioTrack
{
    class FLACDecoder : public FLACAPI::Decoder::File
    {
    public:
        FLACDecoder() : FLACAPI::Decoder::File()
        {
            out = nullptr;
        }

        ~FLACDecoder()
        {
            if (out != nullptr) delete out;
        }

        // return constructed flac. call after read done.
        FLACTrack get();

    private:

        virtual FLAC__StreamDecoderWriteStatus write_callback(
            const FLAC__Frame* frame,
            const FLAC__int32* const* buffer) override;

        virtual void metadata_callback(
            const FLAC__StreamMetadata* metadata) override;

        virtual void error_callback(
            FLAC__StreamDecoderErrorStatus status) override;

        FLACTrack* out;
        bool finished = false;
        unsigned int bit_depth = 0; // of samples
    };

public:

    FLACTrack() : AudioTrack() {}
    FLACTrack(unsigned long length);
    FLACTrack(unsigned long length, unsigned long sample_rate);

    bool readFile(std::string path);
    bool writeFile(std::string path);

};

FLACTrack::FLACTrack(unsigned long length) : AudioTrack(length)
{
    // flac specific info TODO
}

FLACTrack::FLACTrack(unsigned long length, unsigned long sample_rate)
    : AudioTrack(length, sample_rate)
{
}

bool FLACTrack::readFile(std::string path)
{
    FILE* fp = fopen(path.c_str(), "rb");
    if (fp == NULL) throw std::runtime_error("file not found");

    FLACDecoder decoder;
    if (!decoder) throw 0xF;

    decoder.set_md5_checking(true);
    FLAC__StreamDecoderInitStatus init_stat = decoder.init(fp);
    if (init_stat != FLAC__STREAM_DECODER_INIT_STATUS_OK)
        throw std::runtime_error("file read error");

    if (!decoder.process_until_end_of_stream())
        throw std::runtime_error("file read error");

    fclose(fp);

    *this = decoder.get();

    return false;
}

bool FLACTrack::writeFile(std::string path)
{
    return false;
}


FLACTrack FLACTrack::FLACDecoder::get()
{
    if (!finished) throw std::runtime_error("get flac before finished decode");
    return *out;
}

FLAC__StreamDecoderWriteStatus FLACTrack::FLACDecoder::write_callback(const FLAC__Frame* frame, const FLAC__int32* const* buffer)
{
    // ASSUMES out is alr. constructed
    // TODO add variable stuff support
    if (out->sampleLength() == 0) throw std::runtime_error("sample length must be > 0");
    if (out->channels() != 2) throw std::invalid_argument("unsupported");

    //const FLAC__uint32 total_size = (FLAC__uint32)(out->sampleLength() * out->channels() * (bit_depth / 8));

    for (unsigned long i = 0; i < frame->header.blocksize; i++)
    {
        audio_sample s{};
        s.left = (float)buffer[0][i] / uint_max_val(bit_depth); // bitdepth max int
        s.right = (float)buffer[1][i] / uint_max_val(bit_depth);
        out->audio[out->currentPos()] = s;
        (void)out->nextSample();
    }

    finished = true;
    return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE;
}

void FLACTrack::FLACDecoder::metadata_callback(const FLAC__StreamMetadata* metadata)
{
    if (metadata->type != FLAC__METADATA_TYPE_STREAMINFO) return;

    unsigned long sample_len =
        (unsigned long)metadata->data.stream_info.total_samples;
    unsigned long sample_rate =
        (unsigned long)metadata->data.stream_info.sample_rate;
    unsigned int channels = metadata->data.stream_info.channels;

    if (channels != 2) throw std::invalid_argument("unsupported");
    out = new FLACTrack(sample_len, sample_rate);

    bit_depth = metadata->data.stream_info.bits_per_sample;
}

void FLACTrack::FLACDecoder::error_callback(FLAC__StreamDecoderErrorStatus status)
{
    throw std::runtime_error("FLAC decoder error: " + std::to_string(status));
    return;
}

class AudioPlayer
{
public:
    AudioPlayer();
    ~AudioPlayer();

    virtual void operator()(AudioTrack& track);
    void play(AudioTrack& track);
    void stop();
    void close();
    AudioTrack* track();
    void setGain(float db); // in db (10log_10(out/in))

private:

    static int playCallback(
        const void* inputBuffer,
        void* outputBuffer,
        unsigned long framesPerBuffer,
        const PaStreamCallbackTimeInfo* timeInfo,
        PaStreamCallbackFlags statusFlags,
        void* userData);

    bool handle_pa_error(PaError err);
    PaStream* stream; /**< Port Audio stream */
    AudioTrack* audio_track; /**< Audio track to play */
    float amplifier; /**< Linear gain scalar */

};

AudioPlayer::AudioPlayer() : amplifier(1.0)
{
    stream = nullptr;
    PaError err;
    err = Pa_Initialize();
    handle_pa_error(err);
}

AudioPlayer::~AudioPlayer()
{
    close();
    handle_pa_error(Pa_Terminate());
}

void AudioPlayer::operator()(AudioTrack& track)
{
    close();
    play(track);
}

void AudioPlayer::play(AudioTrack& track)
{
    audio_track = &track;
    PaError err;
    /* Open an audio I/O stream. */
    err = Pa_OpenDefaultStream(&stream,
        0,          /* no input channels */
        2,          /* stereo output */
        paFloat32,  /* 32 bit floating point output */
        track.sampleRate(),
        paFramesPerBufferUnspecified, /* frames per buffer, i.e. the number
                           of sample frames that PortAudio will
                           request from the callback. Many apps
                           may want to use
                           paFramesPerBufferUnspecified, which
                           tells PortAudio to pick the best,
                           possibly changing, buffer size.*/
        playCallback, /* this is your callback function */
        this); /*This is a pointer that will be passed to
                           your callback*/
    handle_pa_error(err);

    track.resetPosition();
    err = Pa_StartStream(stream);
    handle_pa_error(err);

    //Pa_Sleep(1000); // handle duration?
    while (Pa_IsStreamActive(stream)) {}
}

void AudioPlayer::stop()
{ /// AbortStream() for immediate stop (forceq)
    if (stream == nullptr) return; // nothing to stop
    handle_pa_error(Pa_StopStream(stream));
}


void AudioPlayer::close()
{
    if (stream == nullptr) return; // nothing to close
    stop();
    PaError err;
    err = Pa_CloseStream(stream);
    handle_pa_error(err);

}

AudioTrack* AudioPlayer::track()
{
    return audio_track;
}

void AudioPlayer::setGain(float db)
{
    amplifier = pow(10.0, db / 10.0); // db to ratio
}

int AudioPlayer::playCallback(
    const void* inputBuffer, // track audio array
    void* outputBuffer,
    unsigned long framesPerBuffer,
    const PaStreamCallbackTimeInfo* timeInfo,
    PaStreamCallbackFlags statusFlags,
    void* userData)
{
    /* Cast data passed through stream to our structure. */
    AudioPlayer* player = (AudioPlayer*)userData;
    if (!player) throw 0xFF; // if actually happends then WTF??
    AudioTrack* track = player->track();
    if (!track) return paAbort; // cant play nothing.
    audio_sample data = track->nextSample();

    float* out = (float*)outputBuffer;
    (void)inputBuffer; /* Prevent unused variable warning. */

    for (unsigned int i = 0; i < framesPerBuffer; i++)
    {
        if (track->ended()) return paComplete;
        *out++ = player->amplifier * (float)data.left;  /* left */
        *out++ = player->amplifier * (float)data.right;  /* right */
        data = track->nextSample();
    }

    return paContinue;
}

bool AudioPlayer::handle_pa_error(PaError err)
{
    if (err != paNoError) {
        Pa_Terminate();
        throw err;
    }
    return true;
}


class RGBAPixel
{
public:
    RGBAPixel();
    RGBAPixel(uint8_t r, uint8_t g, uint8_t b);
    RGBAPixel(uint8_t r, uint8_t g, uint8_t b, uint8_t a);

    void operator=(const RGBAPixel& other);
    bool operator==(const RGBAPixel& other) const;
    bool operator!=(const RGBAPixel& other) const;
    void set(uint8_t r, uint8_t g, uint8_t b);
    void set(uint8_t r, uint8_t g, uint8_t b, uint8_t a);
    void set(const RGBAPixel& other);
    void setAlpha(uint8_t a);
    void setTransparency(float transparency);
    float dist(const RGBAPixel& other) const;
    float intensity() const;

    uint8_t r;
    uint8_t g;
    uint8_t b;
    uint8_t a;

private:

    void clampRGBA();

};


RGBAPixel::RGBAPixel()
    : RGBAPixel(0, 0, 0, 0)
{

}

RGBAPixel::RGBAPixel(uint8_t r, uint8_t g, uint8_t b)
    : RGBAPixel(r, g, b, 255)
{
}

RGBAPixel::RGBAPixel(uint8_t r, uint8_t g, uint8_t b, uint8_t a) : r(r), g(g), b(b), a(a)
{
    clampRGBA();
}

void RGBAPixel::operator=(const RGBAPixel& other)
{
    this->r = other.r;
    this->g = other.g;
    this->b = other.b;
    this->a = other.a;
    clampRGBA();
}

bool RGBAPixel::operator==(const RGBAPixel& other) const
{
    return r == other.r && g == other.g && b == other.b && a == other.a;
}

bool RGBAPixel::operator!=(const RGBAPixel& other) const
{
    return !(*this == other);
}


// set rgb
void RGBAPixel::set(uint8_t r, uint8_t g, uint8_t b)
{
    *this = RGBAPixel(r, g, b, this->a);
}


void RGBAPixel::set(uint8_t r, uint8_t g, uint8_t b, uint8_t a)
{
    *this = RGBAPixel(r, g, b, a);
}

void RGBAPixel::set(const RGBAPixel& other)
{
    *this = other;
}

void RGBAPixel::setAlpha(uint8_t a)
{
    this->a = a;
    clampRGBA();
}

void RGBAPixel::setTransparency(float transparency)
{
    setAlpha(1 - transparency);
}

// return euclidian distance 
float RGBAPixel::dist(const RGBAPixel& other) const
{
    float r_dist = r - other.r;
    float g_dist = g - other.g;
    float b_dist = b - other.b;
    return sqrt(r_dist * r_dist + g_dist * g_dist + b_dist * b_dist);
}

float RGBAPixel::intensity() const
{
    return ((float) r + g + b) / (3.0 * 255);
}

// clamp all [0,1]
// Assuming NaN reprsents infinity, clamp to 1.
void RGBAPixel::clampRGBA()
{
    r = std::max(std::min(r, (uint8_t)255), (uint8_t)0);
    g = std::max(std::min(g, (uint8_t)255), (uint8_t)0);
    b = std::max(std::min(b, (uint8_t)255), (uint8_t)0);
    a = std::max(std::min(a, (uint8_t)255), (uint8_t)0);
}


class PNG
{
    struct pixel_array_struct
    {
        RGBAPixel* pixels;
        int width;
        int height;
    };

public:
    PNG() 
    {
        pixel_data.width = 0;
        pixel_data.height = 0;
        pixel_data.pixels = nullptr;
    }
    ~PNG() { clear(); }

    bool readFile(std::string path);
    void writeFile(std::string path);
    void clear();
    void setInitialSize(int w, int h);
    int width() const;
    int height() const;
    int size() const;
    bool setRGBAPixel(const RGBAPixel& color, int x, int y);
    RGBAPixel* getRGBAPixel(int x, int y) const;
    RGBAPixel* pixelAt(int x, int y) const;

private:
    bool isFilePNG(FILE* fp);

    bool create_png_read_structs(png_structp& png_p, png_infop& info_p, png_infop& end_info_p);
    bool create_png_write_structs(png_structp& png_p, png_infop& info_p);
    int process_png_info(png_structp& png_p, png_infop& info_p);
    int write_png_info(png_structp& png_p, png_infop& info_p);
    int process_png_pixels(png_structp& png_p, png_infop& info_p);
    int write_png_pixels(png_structp& png_p, png_infop& info_p);
    void row_ptrs_to_pixels(const png_bytepp& row_ptrs);
    void pixels_to_row_ptrs(png_bytepp& row_ptrs);

    pixel_array_struct pixel_data;
};

void PNG::clear()
{
    if (pixel_data.pixels != nullptr)
    {
        delete[] pixel_data.pixels;
        pixel_data.pixels = nullptr;
    }
}

void PNG::setInitialSize(int w, int h)
{
    if (w <= 0 || h <= 0) throw std::invalid_argument("image size must be positive integer");
    if (w > 100000 || h > 100000) throw std::invalid_argument("image too large"); // artificial limit
    pixel_data.width = w;
    pixel_data.height = h;
    pixel_data.pixels = new RGBAPixel[size()];
}

int PNG::height() const
{
    return pixel_data.height;
}

int PNG::width() const
{
    return pixel_data.width;
}

int PNG::size() const
{
    return width() * height();
}

bool PNG::setRGBAPixel(const RGBAPixel& color, int x, int y)
{
    RGBAPixel* target = pixelAt(x, y);
    if (target == nullptr) return false;
    *target = color;
    return true;
}

RGBAPixel* PNG::getRGBAPixel(int x, int y) const
{
    RGBAPixel* p = pixelAt(x, y);
    if (p == nullptr) throw std::out_of_range("Image:: Pixel index out of range: (" +
        std::to_string(x) + ", " +
        std::to_string(y) + ")");
    return p;
}

RGBAPixel* PNG::pixelAt(int x, int y) const
{
    if (pixel_data.pixels == NULL) return nullptr;
    if (x < 0 || x >= width() || y < 0 || y >= height()) return nullptr;
    return &pixel_data.pixels[y * width() + x];
}

bool PNG::readFile(std::string path)
{
    clear();
    FILE* fp = fopen(path.c_str(), "rb"); // read-binary file

    try {
        if (!isFilePNG(fp)) throw std::invalid_argument("file " + path + " is not a png");
    }
    catch (const std::exception &e)
    {
        return false;
    }
    // png structs low level interface libpng: Create read, info and end info structs
    png_structp png_ptr;
    png_infop info_ptr, end_info_ptr;
    if (!create_png_read_structs(png_ptr, info_ptr, end_info_ptr))
    {
        throw std::runtime_error("Failed to read png");
    }
    png_set_sig_bytes(png_ptr, PNG_HEADER_BYTES);
    // Set false jump buffer
    if (setjmp(png_jmpbuf(png_ptr)))
    {
        png_destroy_read_struct(&png_ptr, &info_ptr, &end_info_ptr);
        if (fp) fclose(fp);
        throw std::runtime_error("Failed to read png");
    }
    png_init_io(png_ptr, fp);
    // read info to class PNG
    if (!process_png_info(png_ptr, info_ptr)) throw std::runtime_error("Failed to read png");
    if (!process_png_pixels(png_ptr, info_ptr)) throw std::runtime_error("Failed to read png");
    png_destroy_read_struct(&png_ptr, &info_ptr, &end_info_ptr);
    if (fp) fclose(fp);

    return true;
}

void PNG::writeFile(std::string path)
{
    FILE* fp = fopen(path.c_str(), "wb");
    if (fp == NULL)
    {
        throw std::runtime_error("Failed to write png");
    }

    png_structp write_png_ptr;
    png_infop write_info_ptr;
    if (!create_png_write_structs(write_png_ptr, write_info_ptr)) throw std::runtime_error("Failed to write png");


    png_init_io(write_png_ptr, fp);
    if (!write_png_info(write_png_ptr, write_info_ptr)) throw std::runtime_error("Failed to write png");
    if (!write_png_pixels(write_png_ptr, write_info_ptr)) throw std::runtime_error("Failed to write png");

    png_destroy_write_struct(&write_png_ptr, &write_info_ptr);
    fclose(fp);
}

// return true if file is PNG using libpng
bool PNG::isFilePNG(FILE* fp)
{
    if (!fp) throw std::invalid_argument("File does not exist");
    png_byte header[PNG_HEADER_BYTES];
    fread(header, 1, PNG_HEADER_BYTES, fp);
    bool is_png = !png_sig_cmp(header, 0, PNG_HEADER_BYTES);
    if (!is_png) return false;
    return true;
}

// creates png and info structs from libpng
bool PNG::create_png_read_structs(png_structp& png_p, png_infop& info_p, png_infop& end_info_p)
{
    // could add false and warning handlers here
    png_p = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_p) return false;
    info_p = png_create_info_struct(png_p);
    end_info_p = png_create_info_struct(png_p);
    if (!info_p || !end_info_p)
    {
        png_destroy_read_struct(&png_p, &info_p, &end_info_p);
        return false;
    }
    return true;
}

bool PNG::create_png_write_structs(png_structp& png_p, png_infop& info_p) // optional: end_write_info_ptr;
{
    png_p = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_p) return false;
    info_p = png_create_info_struct(png_p);
    if (!info_p)
    {
        png_destroy_write_struct(&png_p, NULL);
        return false;
    }
    return true;
}

// process info into PNG class (this)
int PNG::process_png_info(png_structp& png_p, png_infop& info_p)
{
    png_read_info(png_p, info_p);
    png_uint_32 width, height;
    int bit_depth, color_type, filter, compression, interlace;
    png_get_IHDR(png_p, info_p, &width, &height, &bit_depth, &color_type, &filter, &compression, &interlace);
    setInitialSize(width, height);
    png_set_expand(png_p); // expand to rgba8, same as below:
    png_set_add_alpha(png_p, 0xFFFF, PNG_FILLER_AFTER);
    png_set_strip_16(png_p);

    png_read_update_info(png_p, info_p);
    return true;
}

// process metadata, compression etc, from PNG class (this)
int PNG::write_png_info(png_structp& png_p, png_infop& info_p)
{
    png_set_IHDR(png_p, info_p, width(), height(), 8,
        PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_set_filter(png_p, 0, PNG_FILTER_NONE);
    png_set_compression_level(png_p, PNG_Z_DEFAULT_COMPRESSION);
    png_time p_time;
    time_t ttime;
    time(&ttime);
    png_convert_from_time_t(&p_time, ttime);
    png_set_tIME(png_p, info_p, &p_time);
    //png_set_text(write_png_p, write_info_p, <some png_text>, <how many texts>);
    png_write_info(png_p, info_p);

    return true;
}

// read png_bytes to PNG, handle memory operations
int PNG::process_png_pixels(png_structp& png_p, png_infop& info_p)
{
    int status = true;
    // read png bytes (img data) - alloc memory for libpng copy then transfer to PNG class (this)
    png_byte channels = png_get_channels(png_p, info_p);
    png_byte bit_depth = png_get_bit_depth(png_p, info_p) / 8;
    // alloc height * space per row ptr
    png_bytep* row_ptrs = (png_bytep*)png_malloc(png_p, height() * sizeof(png_bytep));
    for (unsigned int y = 0; y < height(); y++)
    {
        // alloc width * space per pixel // TODO bit depth != 8
        row_ptrs[y] = (png_byte*)png_malloc(png_p, width() * sizeof(png_byte) * channels);
    }
    png_read_image(png_p, row_ptrs);
    // transfer row ptrs to PNG class
    row_ptrs_to_pixels(row_ptrs);
    // free row ptrs after finished
    for (unsigned int y = 0; y < height(); y++)
    {
        png_free(png_p, row_ptrs[y]);
    }
    png_free(png_p, row_ptrs);

    return status;
}

// write PNG bitmap to png_bytes
int PNG::write_png_pixels(png_structp& png_p, png_infop& info_p)
{
    int status = true;

    png_bytepp row_ptrs = (png_bytep*)png_malloc(png_p, height() * sizeof(png_bytep));
    for (unsigned int y = 0; y < height(); y++)
    {
        row_ptrs[y] = (png_byte*)png_malloc(png_p, width() * sizeof(png_byte) * 4); // for not don't handle bit depth != 8
    }
    pixels_to_row_ptrs(row_ptrs);
    png_set_rows(png_p, info_p, row_ptrs);

    png_write_png(png_p, info_p, PNG_TRANSFORM_IDENTITY, NULL); 
    for (unsigned int y = 0; y < height(); y++)
    {
        png_free(png_p, row_ptrs[y]);
    }
    png_free(png_p, row_ptrs);
    row_ptrs = NULL;


    return status;
}

void PNG::row_ptrs_to_pixels(const png_bytepp& row_ptrs)
{
    bool alpha = true; png_byte channels = 4;
    for (unsigned int y = 0; y < height(); y++)
    {
        for (unsigned int x = 0; x < width(); x++)
        {
            unsigned int pixel_size = channels * sizeof(png_byte);
            png_uint_32 r = 0, g = 0, b = 0, a = 0;
            png_uint_32 base = x * pixel_size;
            png_bytep byte_data = &row_ptrs[y][base];
            r = byte_data[0];
            g = byte_data[1];
            b = byte_data[2];
            a = byte_data[3];
            RGBAPixel p(r, g, b, a); 
            setRGBAPixel(p, x, y);
        }
    }
}

void PNG::pixels_to_row_ptrs(png_bytepp& row_ptrs)
{
    // up to 16 bit per channel  TODO, use PNG class fields
    bool alpha = true; png_byte bit_depth = 8;
    png_byte channels = alpha ? 4 : 3; // force alpha for now
    for (unsigned int y = 0; y < height(); y++)
    {
        for (unsigned int x = 0; x < width(); x++)
        {
            unsigned int pixel_size = channels * sizeof(png_byte);
            png_uint_32 base = x * pixel_size;
            png_bytep byte_data = &row_ptrs[y][base];
            // convert all to int
            RGBAPixel p = *getRGBAPixel(x, y);
            png_byte r, g, b, a;
            r = (png_byte)p.r;
            g = (png_byte)p.g;
            b = (png_byte)p.b;
            a = (png_byte)p.a;
            byte_data[0] = r;
            byte_data[1] = g;
            byte_data[2] = b;
            byte_data[3] = a;
        }
    }
}

class Timer
{
public:
    Timer();

    void reset();
    void start();

    template <class T = double, class Period = std::ratio<1>>
    T lap()
    {
        using namespace std::chrono;
        duration<T, Period> elapsed;
        if (running) {
            high_resolution_clock::time_point now_point = high_resolution_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::duration<T, Period>>(now_point - last_lap_point);
            last_lap_point = now_point;
        }
        else {
            elapsed = std::chrono::duration_cast<std::chrono::duration<T, Period>>(stop_point - last_lap_point);
        }
        return elapsed.count();
    }

    template <class T = double, class Period = std::ratio<1>>
    T time()
    {
        using namespace std::chrono;
        duration<T, Period> elapsed;
        if (running) {
            high_resolution_clock::time_point now_point = high_resolution_clock::now();
            elapsed = std::chrono::duration_cast<std::chrono::duration<T, Period>>(now_point - start_point);
        }
        else {
            elapsed = std::chrono::duration_cast<std::chrono::duration<T, Period>>(stop_point - start_point);
        }
        return elapsed.count();
    }

    template <class T = double, class Period = std::ratio<1>>
    T stop()
    {
        using namespace std::chrono;
        if (running) stop_point = high_resolution_clock::now();
        running = false;
        return time<T, Period>();
    }

private:
    std::chrono::high_resolution_clock::time_point start_point;
    std::chrono::high_resolution_clock::time_point last_lap_point;
    std::chrono::high_resolution_clock::time_point stop_point;
    bool running;

};

Timer::Timer()
    : running(false)
{
    reset();
}

void Timer::reset()
{
    using namespace std::chrono;

    last_lap_point = high_resolution_clock::now();
    start_point = last_lap_point;
    stop_point = start_point;
}

void Timer::start()
{
    reset();
    running = true;
}

std::string get_n_img_path(size_t n)
{
    static constexpr const char* img_dir = "../img/";
    char img_num_str[7];
    snprintf(img_num_str, 7, "%06zu", n);
    std::string img_path;
    // _XXXXXX.png
    img_path.reserve(7 + 1 + 6 + 4);
    img_path += img_dir;
    img_path += "_";
    img_path += img_num_str;
    img_path += ".png";

    return img_path;
}

#include <windows.h>

CONSOLE_SCREEN_BUFFER_INFO csbi;
const BOOL gcsbi_stat = GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
const int term_cols = csbi.srWindow.Right - csbi.srWindow.Left + 1;
const int term_rows = csbi.srWindow.Bottom - csbi.srWindow.Top + 1;

const size_t term_buf_size = (term_cols + 1) * term_rows + 1;
char* term_buf1 = new char[term_buf_size];
char*& active_buf = term_buf1;

#include <windows.h>

#ifndef ENABLE_VIRTUAL_TERMINAL_PROCESSING
#define ENABLE_VIRTUAL_TERMINAL_PROCESSING 0x0004
#endif // ENABLE_VIRTUAL_TERMINAL_PROCESSING

#define ENABLE_ANSI() \
        do { \
            HANDLE hout = GetStdHandle(STD_OUTPUT_HANDLE); \
            DWORD dw; \
            GetConsoleMode(hout, &dw); \
            SetConsoleMode(hout, dw | ENABLE_VIRTUAL_TERMINAL_PROCESSING); \
        } while (0)

#define DISABLE_ANSI() \
        do { \
            HANDLE hout = GetStdHandle(STD_OUTPUT_HANDLE); \
            DWORD dw; \
            GetConsoleMode(hout, &dw); \
            SetConsoleMode(hout, dw & ~ENABLE_VIRTUAL_TERMINAL_PROCESSING); \
        } while(0)

int setup()
{
    assert(term_cols > 1);
    assert(term_rows > 1);
    ENABLE_ANSI();
    for (int i = 0; i < term_rows; i++)
    {
        term_buf1[term_buf_size - 1] = '\0';
    }
    return 0;
}

inline char* get_term_row_ptr(char* term, int y)
{
    return &term[y * (term_cols + 1)];
}

inline char* get_active_row_ptr(int y)
{
    return get_term_row_ptr(active_buf, y);
}

inline char* get_active_buf()
{
    return active_buf;
}

char select_char(const RGBAPixel& p)
{
    const float intensity = p.intensity();
    char c = 0;
    if (intensity < 0.1)
        c = ' ';
    else if (intensity < 0.2)
        c = '.';
    else if (intensity < 0.4)
        c = '+';
    else if (intensity < 0.6)
        c = '*';
    else if (intensity < 0.8)
        c = '&';
    else if (intensity < 0.9)
        c = '#';
    else
        c = '@';
    return c;
}

bool check_file_readable(std::string path)
{
    FILE* f = fopen(path.c_str(), "r");
    bool rv = (f != nullptr);
    if (rv) fclose(f);
    return rv;
}

// image load buffer
PNG img;

static int dummy_setup = setup();

PNG* get_png(size_t i)
{
    if (!img.readFile(get_n_img_path(i))) return nullptr;
    return &img;
}

void clear()
{
    std::cout << "\x1B[H" << std::flush;
}

void play_realtime()
{
    PNG* pimg = nullptr;
    size_t i = 0; // which img number to use
    Timer t;
    t.start();
    constexpr size_t target_frametime = std::micro::den / 30; // n microsec
    size_t frame_time = 0;
    size_t time_acc = 0;
    // buffer for sprintf so i can ignore '\0'
    char* last_row_tmp = new char[term_cols];
    // main loop
    while (pimg = get_png(i))
    {
        frame_time = t.lap<size_t, std::micro>();
        time_acc += frame_time;
        PNG& img = *pimg;
        clear();
        for (int y = 0; y < term_rows; y++)
        {
            char* row = get_active_row_ptr(y);
            for (int x = 0; x < term_cols; x++)
            {
                const int pick_x = (int)((double)x * img.width() / term_cols);
                const int pick_y = (int)((double)y * img.height() / term_rows);
                RGBAPixel* p = img.getRGBAPixel(pick_x, pick_y);
                char c = select_char(*p);
                row[x] = c;
            }
            row[term_cols] = '\n';
        }
        // show frame info on last line
        char* last_row = get_active_row_ptr(term_rows - 1);
        // remove \n on last row
        last_row[term_cols] = '\0';
        // printf frame info into bottom left
        size_t info_slen = snprintf(last_row_tmp, term_cols, "%zu %zums %ffps %fs", i, frame_time/1000, (double)i / t.time(), t.time());
        memcpy(last_row, last_row_tmp, min(info_slen, term_cols));
        // show
        std::cout << get_active_buf();
        // keep at 30 fps
        if (frame_time < target_frametime)
        {
            std::this_thread::sleep_for(std::chrono::duration<size_t, std::micro>((target_frametime - frame_time)));
        }
        size_t di = time_acc / target_frametime;
        i += di;
        time_acc = time_acc % target_frametime;
    }
    delete[] last_row_tmp;
}

int main()
{
    AudioPlayer player;
    FLACTrack audio;
    audio.readFile("../audio/audio.flac");
    std::thread audio_thread([&]() -> void
        {
            player.play(audio);
        });
    std::thread render_thread([]() -> void
        {
            play_realtime();
        });
    audio_thread.join();
    render_thread.join();

    return 0;
}