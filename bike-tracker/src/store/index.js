import Vue from 'vue'
import Vuex from 'vuex'

Vue.use(Vuex)

export default new Vuex.Store({
  state: {
    stations: [],
    showAppInfo: true,
    showStationInfo: false,
    onlyAvailable: false,
    showFull: false
  },
  mutations: {
    setStations: (state, stations) => {
      stations.forEach(element => {
        state.stations.push({
          station_id: element.station_id,
          name: element.name,
          pos: {
            lat: element.lat,
            lng: element.lon,
          },
          capacity: element.capacity
        })
      });
    },
    toggleAppInfo: (state) => state.showAppInfo = !state.showAppInfo,
    toggleStationInfo: (state) => state.showStationInfo = !state.showStationInfo,
    toggleDisplay: (state) => state.onlyAvailable = !state.onlyAvailable,
    toggleFull: (state) => state.showFull = !state.showFull
  },    
  actions: {
    async getStations({ commit }) {
      await fetch('https://gbfs.urbansharing.com/oslobysykkel.no/station_information.json').then(
        res => res.json()
      ).then(
         data => {
           commit('setStations', data.data.stations)
           //console.log(data)
           //console.log(this.state.stations)
          }).catch(error => console.log(error))
      //curl -H "Client-Identifier: privat-Citybiketracker" \ https://gbfs.urbansharing.com/oslobysykkel.no/station_information.json
    },
    async getStationData() {
      await fetch('https://gbfs.urbansharing.com/oslobysykkel.no/station_status.json').then(
        res => res.json()
      ).then( data => {
        data.data.stations.forEach(element => {
          var station = this.state.stations.find(station => {return station.station_id === element.station_id})
          station.available = element.num_bikes_available
          station.docksAvailable = element.num_docks_available
        })
      })
    },
    toggleAppInfo({ commit }) {
      commit('toggleAppInfo')
    },
    toggleStationInfo({ commit }) {
      commit('toggleStationInfo')
    },
    toggleDisplay({ commit }) {
      commit('toggleDisplay')
    },
    toggleFull({ commit }) {
      commit('toggleFull')
    }
  },
  modules: {
  },
  getters: {
    stationsFiltered: (state) => {
      //fiks filter
      var returnStat
      if (!state.onlyAvailable && !state.showFull) {
        returnStat = state.stations.filter(stat => { return stat.available > 0 })
      } else if (state.onlyAvailable && !state.showfull) {
        returnStat = state.stations
        //returnStat = state.stations.filter(stat => { return stat.available < stat.capacity })
      } else if (!state.onlyAvailable && state.showFull){
        returnStat = state.stations.filter(stat => { return stat.available < stat.capacity && stat.available > 0 })
      } else if (state.onlyAvailable && state.showFull) {
        returnStat = state.stations.filter(stat => { return stat.available < stat.capacity })
      } else {
        returnStat = state.stations
      }
      returnStat.forEach(stat => console.log(stat.name + ": " + (stat.available - stat.capacity)))
      return returnStat
    },
    showAppInfo: (state) => state.showAppInfo,
    showStationInfo: (state) => state.showStationInfo,
    showUnavailable: (state) => state.onlyAvailable,
    showFull: (state) => state.showFull
  }
})