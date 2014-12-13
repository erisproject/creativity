--
-- PostgreSQL database dump
--

SET statement_timeout = 0;
SET lock_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SET check_function_bodies = false;
SET client_min_messages = warning;

SET search_path = public, pg_catalog;

ALTER TABLE IF EXISTS ONLY public.state DROP CONSTRAINT IF EXISTS state_sim_fkey;
ALTER TABLE IF EXISTS ONLY public.setting DROP CONSTRAINT IF EXISTS setting_sim_fkey;
ALTER TABLE IF EXISTS ONLY public.reader DROP CONSTRAINT IF EXISTS reader_state_fkey;
ALTER TABLE IF EXISTS ONLY public.library DROP CONSTRAINT IF EXISTS library_reader_fkey;
ALTER TABLE IF EXISTS ONLY public.library DROP CONSTRAINT IF EXISTS library_book_fkey;
ALTER TABLE IF EXISTS ONLY public.friend DROP CONSTRAINT IF EXISTS friend_reader_fkey;
ALTER TABLE IF EXISTS ONLY public.friend DROP CONSTRAINT IF EXISTS friend_friend_fkey;
ALTER TABLE IF EXISTS ONLY public.book DROP CONSTRAINT IF EXISTS book_state_fkey;
ALTER TABLE IF EXISTS ONLY public.book DROP CONSTRAINT IF EXISTS book_author_fkey;
ALTER TABLE IF EXISTS ONLY public.belief DROP CONSTRAINT IF EXISTS belief_reader_fkey;
ALTER TABLE IF EXISTS ONLY public.state DROP CONSTRAINT IF EXISTS state_sim_t_key;
ALTER TABLE IF EXISTS ONLY public.state DROP CONSTRAINT IF EXISTS state_pkey;
ALTER TABLE IF EXISTS ONLY public.simulation DROP CONSTRAINT IF EXISTS simulation_seed_key;
ALTER TABLE IF EXISTS ONLY public.simulation DROP CONSTRAINT IF EXISTS simulation_pkey;
ALTER TABLE IF EXISTS ONLY public.setting DROP CONSTRAINT IF EXISTS setting_pkey;
ALTER TABLE IF EXISTS ONLY public.reader DROP CONSTRAINT IF EXISTS reader_state_eris_id_key;
ALTER TABLE IF EXISTS ONLY public.reader DROP CONSTRAINT IF EXISTS reader_pkey;
ALTER TABLE IF EXISTS ONLY public.library DROP CONSTRAINT IF EXISTS library_pkey;
ALTER TABLE IF EXISTS ONLY public.friend DROP CONSTRAINT IF EXISTS friend_pkey;
ALTER TABLE IF EXISTS ONLY public.book DROP CONSTRAINT IF EXISTS book_state_eris_id_key;
ALTER TABLE IF EXISTS ONLY public.book DROP CONSTRAINT IF EXISTS book_pkey;
ALTER TABLE IF EXISTS ONLY public.belief DROP CONSTRAINT IF EXISTS belief_pkey;
ALTER TABLE IF EXISTS public.state ALTER COLUMN id DROP DEFAULT;
ALTER TABLE IF EXISTS public.simulation ALTER COLUMN id DROP DEFAULT;
ALTER TABLE IF EXISTS public.reader ALTER COLUMN id DROP DEFAULT;
ALTER TABLE IF EXISTS public.book ALTER COLUMN id DROP DEFAULT;
DROP SEQUENCE IF EXISTS public.state_id_seq;
DROP TABLE IF EXISTS public.state;
DROP SEQUENCE IF EXISTS public.simulation_id_seq;
DROP TABLE IF EXISTS public.simulation;
DROP TABLE IF EXISTS public.setting;
DROP SEQUENCE IF EXISTS public.reader_id_seq;
DROP TABLE IF EXISTS public.reader;
DROP TABLE IF EXISTS public.library;
DROP TABLE IF EXISTS public.friend;
DROP SEQUENCE IF EXISTS public.book_id_seq;
DROP TABLE IF EXISTS public.book;
DROP TABLE IF EXISTS public.belief;
DROP TYPE IF EXISTS public.library_type;
DROP TYPE IF EXISTS public.belief_type;
DROP EXTENSION IF EXISTS plpgsql;
DROP SCHEMA IF EXISTS public;
--
-- Name: public; Type: SCHEMA; Schema: -; Owner: -
--

CREATE SCHEMA public;


--
-- Name: SCHEMA public; Type: COMMENT; Schema: -; Owner: -
--

COMMENT ON SCHEMA public IS 'standard public schema';


--
-- Name: plpgsql; Type: EXTENSION; Schema: -; Owner: -
--

CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;


--
-- Name: EXTENSION plpgsql; Type: COMMENT; Schema: -; Owner: -
--

COMMENT ON EXTENSION plpgsql IS 'PL/pgSQL procedural language';


SET search_path = public, pg_catalog;

--
-- Name: belief_type; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE belief_type AS ENUM (
    'profit',
    'profit_extrap',
    'demand',
    'quality',
    'profit_stream'
);


--
-- Name: library_type; Type: TYPE; Schema: public; Owner: -
--

CREATE TYPE library_type AS ENUM (
    'wrote',
    'bought',
    'pirated'
);


SET default_tablespace = '';

SET default_with_oids = false;

--
-- Name: belief; Type: TABLE; Schema: public; Owner: -; Tablespace: 
--

CREATE TABLE belief (
    reader bigint NOT NULL,
    type belief_type NOT NULL,
    k smallint NOT NULL,
    noninformative boolean NOT NULL,
    s2 double precision,
    n double precision,
    beta double precision[],
    v_lower double precision[]
);


--
-- Name: book; Type: TABLE; Schema: public; Owner: -; Tablespace: 
--

CREATE TABLE book (
    id bigint NOT NULL,
    state integer NOT NULL,
    eris_id bigint NOT NULL,
    author_eris_id bigint NOT NULL,
    created integer NOT NULL,
    "position" double precision[] NOT NULL,
    quality double precision NOT NULL,
    price double precision,
    revenue double precision NOT NULL,
    revenue_lifetime double precision NOT NULL,
    sales integer NOT NULL,
    sales_lifetime integer NOT NULL,
    pirated integer NOT NULL,
    pirated_lifetime integer NOT NULL,
    lifetime integer NOT NULL
);


--
-- Name: book_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE book_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: book_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE book_id_seq OWNED BY book.id;


--
-- Name: friend; Type: TABLE; Schema: public; Owner: -; Tablespace: 
--

CREATE TABLE friend (
    reader bigint NOT NULL,
    friend_eris_id bigint NOT NULL
);


--
-- Name: library; Type: TABLE; Schema: public; Owner: -; Tablespace: 
--

CREATE TABLE library (
    reader bigint NOT NULL,
    book_eris_id bigint NOT NULL,
    type library_type NOT NULL,
    new boolean NOT NULL,
    quality double precision NOT NULL
);


--
-- Name: reader; Type: TABLE; Schema: public; Owner: -; Tablespace: 
--

CREATE TABLE reader (
    id bigint NOT NULL,
    state integer NOT NULL,
    eris_id bigint NOT NULL,
    "position" double precision[] NOT NULL,
    u double precision NOT NULL,
    u_lifetime double precision NOT NULL,
    cost_fixed double precision NOT NULL,
    cost_unit double precision NOT NULL,
    cost_piracy double precision NOT NULL,
    income double precision NOT NULL
);


--
-- Name: reader_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE reader_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: reader_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE reader_id_seq OWNED BY reader.id;


--
-- Name: setting; Type: TABLE; Schema: public; Owner: -; Tablespace: 
--

CREATE TABLE setting (
    sim integer NOT NULL,
    name character varying(100) NOT NULL,
    dbl double precision,
    i64 bigint
);


--
-- Name: simulation; Type: TABLE; Schema: public; Owner: -; Tablespace: 
--

CREATE TABLE simulation (
    id integer NOT NULL,
    seed bigint NOT NULL,
    created timestamp with time zone DEFAULT NOW()
);


--
-- Name: simulation_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE simulation_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: simulation_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE simulation_id_seq OWNED BY simulation.id;


--
-- Name: state; Type: TABLE; Schema: public; Owner: -; Tablespace: 
--

CREATE TABLE state (
    id integer NOT NULL,
    sim integer NOT NULL,
    t integer NOT NULL
);


--
-- Name: state_id_seq; Type: SEQUENCE; Schema: public; Owner: -
--

CREATE SEQUENCE state_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


--
-- Name: state_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: -
--

ALTER SEQUENCE state_id_seq OWNED BY state.id;


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY book ALTER COLUMN id SET DEFAULT nextval('book_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY reader ALTER COLUMN id SET DEFAULT nextval('reader_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY simulation ALTER COLUMN id SET DEFAULT nextval('simulation_id_seq'::regclass);


--
-- Name: id; Type: DEFAULT; Schema: public; Owner: -
--

ALTER TABLE ONLY state ALTER COLUMN id SET DEFAULT nextval('state_id_seq'::regclass);


--
-- Name: belief_pkey; Type: CONSTRAINT; Schema: public; Owner: -; Tablespace: 
--

ALTER TABLE ONLY belief
    ADD CONSTRAINT belief_pkey PRIMARY KEY (reader, type, k);


--
-- Name: book_pkey; Type: CONSTRAINT; Schema: public; Owner: -; Tablespace: 
--

ALTER TABLE ONLY book
    ADD CONSTRAINT book_pkey PRIMARY KEY (id);


--
-- Name: book_state_eris_id_key; Type: CONSTRAINT; Schema: public; Owner: -; Tablespace: 
--

ALTER TABLE ONLY book
    ADD CONSTRAINT book_state_eris_id_key UNIQUE (state, eris_id);


--
-- Name: friend_pkey; Type: CONSTRAINT; Schema: public; Owner: -; Tablespace: 
--

ALTER TABLE ONLY friend
    ADD CONSTRAINT friend_pkey PRIMARY KEY (reader, friend_eris_id);


--
-- Name: library_pkey; Type: CONSTRAINT; Schema: public; Owner: -; Tablespace: 
--

ALTER TABLE ONLY library
    ADD CONSTRAINT library_pkey PRIMARY KEY (reader, book_eris_id);


--
-- Name: reader_pkey; Type: CONSTRAINT; Schema: public; Owner: -; Tablespace: 
--

ALTER TABLE ONLY reader
    ADD CONSTRAINT reader_pkey PRIMARY KEY (id);


--
-- Name: reader_state_eris_id_key; Type: CONSTRAINT; Schema: public; Owner: -; Tablespace: 
--

ALTER TABLE ONLY reader
    ADD CONSTRAINT reader_state_eris_id_key UNIQUE (state, eris_id);


--
-- Name: setting_pkey; Type: CONSTRAINT; Schema: public; Owner: -; Tablespace: 
--

ALTER TABLE ONLY setting
    ADD CONSTRAINT setting_pkey PRIMARY KEY (sim, name);


--
-- Name: simulation_pkey; Type: CONSTRAINT; Schema: public; Owner: -; Tablespace: 
--

ALTER TABLE ONLY simulation
    ADD CONSTRAINT simulation_pkey PRIMARY KEY (id);


--
-- Name: simulation_seed_key; Type: CONSTRAINT; Schema: public; Owner: -; Tablespace: 
--

ALTER TABLE ONLY simulation
    ADD CONSTRAINT simulation_seed_key UNIQUE (seed);


--
-- Name: state_pkey; Type: CONSTRAINT; Schema: public; Owner: -; Tablespace: 
--

ALTER TABLE ONLY state
    ADD CONSTRAINT state_pkey PRIMARY KEY (id);


--
-- Name: state_sim_t_key; Type: CONSTRAINT; Schema: public; Owner: -; Tablespace: 
--

ALTER TABLE ONLY state
    ADD CONSTRAINT state_sim_t_key UNIQUE (sim, t);


--
-- Name: belief_reader_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY belief
    ADD CONSTRAINT belief_reader_fkey FOREIGN KEY (reader) REFERENCES reader(id) ON DELETE CASCADE;


--
-- Name: book_state_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY book
    ADD CONSTRAINT book_state_fkey FOREIGN KEY (state) REFERENCES state(id) ON DELETE CASCADE;


--
-- Name: friend_reader_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY friend
    ADD CONSTRAINT friend_reader_fkey FOREIGN KEY (reader) REFERENCES reader(id) ON DELETE CASCADE;


--
-- Name: library_reader_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY library
    ADD CONSTRAINT library_reader_fkey FOREIGN KEY (reader) REFERENCES reader(id) ON DELETE CASCADE;


--
-- Name: reader_state_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY reader
    ADD CONSTRAINT reader_state_fkey FOREIGN KEY (state) REFERENCES state(id) ON DELETE CASCADE;


--
-- Name: setting_sim_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY setting
    ADD CONSTRAINT setting_sim_fkey FOREIGN KEY (sim) REFERENCES simulation(id) ON DELETE CASCADE;


--
-- Name: state_sim_fkey; Type: FK CONSTRAINT; Schema: public; Owner: -
--

ALTER TABLE ONLY state
    ADD CONSTRAINT state_sim_fkey FOREIGN KEY (sim) REFERENCES simulation(id) ON DELETE CASCADE;


--
-- PostgreSQL database dump complete
--

