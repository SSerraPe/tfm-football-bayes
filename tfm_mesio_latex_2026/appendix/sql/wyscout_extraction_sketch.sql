-- ============================================================================
-- Condensed sketch of the LaLiga player-season extraction and construction
-- query. Reproduces the logic of the seven stages below without the full
-- ~90-column select lists of the actual query, which is kept in this
-- directory (wyscout_extraction_query.sql) for anyone who needs the literal
-- SQL. Grain: one row per player_id x season_id x competition_id.
-- ============================================================================

with

-- 1) Competition and season scope ---------------------------------------------
target_competitions as (
    select competition_id, competition_name, competition_area_name
    from stg_wyscout__competitions
    where lower(competition_name) like '%laliga%'
      and lower(coalesce(competition_gender, 'male')) = 'male'
      and coalesce(competition_division_level, 1) = 1   -- top flight only
),
target_seasons as (
    select season_id, season_name, competition_id
    from stg_wyscout__seasons
    where competition_id in (select competition_id from target_competitions)
),

-- 2) Match-player rows, restricted to the target competition/seasons ---------
match_player_base as (
    select
        mps.player_id, mps.match_id, mps.season_id, mps.competition_id,
        p.player_name, p.role_name,
        mps.total_minutes_on_field,
        -- ~90 raw per-match counting statistics, grouped into blocks:
        --   attacking output      (goals, shots, xG, touches in box, ...)
        --   passing and creation  (passes, key passes, crosses, ...)
        --   carrying / 1v1        (dribbles, progressive runs, ...)
        --   defence and pressing  (interceptions, recoveries, ...)
        --   duels                 (defensive/offensive/aerial duels, won/lost)
        --   discipline            (fouls, cards, losses)
        -- plus a block of goalkeeper-specific markers kept only to drive the
        -- exclusion step below (total_gk_saves, total_gk_shots_against,
        -- total_gk_clean_sheets, total_gk_exits, total_gk_conceded_goals)
        mps.*
    from stg_wyscout__matches_players_stats as mps
    join target_seasons  using (season_id)
    join target_competitions using (competition_id)
    left join stg_wyscout__players as p using (player_id)
),

-- 3) Goalkeeper exclusion, by behaviour rather than by role label -------------
--    Role labels alone are not reliable (outfield players who briefly go in
--    goal, mislabelled positions); goalkeeper-specific activity is not.
outfield_match_player_base as (
    select *
    from match_player_base
    where coalesce(total_gk_shots_against, 0)  = 0
      and coalesce(total_gk_saves, 0)          = 0
      and coalesce(total_gk_conceded_goals, 0) = 0
      and coalesce(total_gk_clean_sheets, 0)   = 0
      and coalesce(total_gk_exits, 0)          = 0
),

-- 4) Aggregate match rows to one row per player x season ----------------------
player_season_base as (
    select
        player_id, season_id, competition_id,
        max(player_name) as player_name,
        sum(total_minutes_on_field) as minutes,
        -- one sum() per raw counting statistic carried over from step 2
        -- (goals, assists, passes, duels, cards, ... ~90 columns)
        sum(total_goals) as goals,
        sum(total_passes) as passes,
        sum(total_duels)  as duels
        -- ... (remaining raw counts omitted here, summed identically)
    from outfield_match_player_base
    group by player_id, season_id, competition_id
),

-- 5) Minimum-exposure filter ---------------------------------------------------
--    Drops small-sample player-seasons whose per-90 rates would otherwise be
--    dominated by noise.
player_season_filtered as (
    select *
    from player_season_base
    where minutes >= 600
),

-- 6) Feature construction: per-90 rates and efficiency ratios ------------------
player_season_features as (
    select
        psf.*,
        -- Per-90 intensity: one term per volume feature (~30 features total),
        -- covering attacking output, passing, carrying, defence, and duels.
        90.0 * psf.goals  / nullif(psf.minutes, 0) as goals_p90,
        90.0 * psf.passes / nullif(psf.minutes, 0) as passes_p90,
        90.0 * psf.duels  / nullif(psf.minutes, 0) as duels_p90,
        -- ... (remaining per-90 features constructed identically)

        -- Efficiency ratios: successful / attempted, one term per ratio
        -- feature (~18 features total: pass accuracy, dribble success rate,
        -- duel win rates by type, shot conversion, goal-per-shot, ...).
        1.0 * psf.successful_passes / nullif(psf.passes, 0) as pass_accuracy,
        1.0 * psf.duels_won         / nullif(psf.duels, 0)  as duel_win_rate
        -- ... (remaining ratio features constructed identically)
    from player_season_filtered as psf
)

-- 7) Final analytical output: one row per player-season -----------------------
select * from player_season_features;
