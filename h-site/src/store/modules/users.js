import db from "../../firebaseConfig.js";

const state = {
  authenticated: false,
  user: {},
};

const getters = {
  auth: (state) => state.authenticated,
  user: (state) => state.user,
};

const actions = {
  async authentication({ commit }, logInfo) {
    await db.auth
      .signInWithEmailAndPassword(logInfo.email, logInfo.password)
      .then(function(user) {
        console.log("Signed in with: ", logInfo.email);
        commit("authenticate", true);
        commit("setUser", user);
      })
      .catch(function(error) {
        console.error("Error signing in: ", error);
      });
  },
  async signUp({ commit }, logInfo) {
    await db.auth
      .createUserWithEmailAndPassword(logInfo.email, logInfo.password)
      .then(function(user) {
        console.log("Signed up with: ", logInfo.email);
        commit("authenticate", true);
        commit("setUser", user);
      })
      .catch(function(error) {
        console.error("Error while signing up: ", error);
      });
  },
  async logOut({ commit }) {
    await db.auth
      .signOut()
      .then(function() {
        console.log("Signout successfull");
        commit("authenticate", false);
      })
      .catch(function(error) {
        console.error("An error occured while signing out: ", error);
      });
  },
  async updateDisplayName(newName) {
    await db.auth.currentUser
      .updateProfile({
        displayName: newName,
      })
      .then(function() {
        console.log("Username set");
      })
      .catch(function(error) {
        console.error("An error accured updating displayName: ", error);
      });
  },
  async deleteAccount({ commit }) {
    await db.auth.currentUser
      .delete()
      .then(function() {
        commit("authenticate", false);
      })
      .catch(function(error) {
        console.error("An error occured deleting account: ", error);
      });
  },
};

const mutations = {
  authenticate: (state, authed) => (state.authenticated = authed),
  setUser: (state, { user }) => (state.user = user),
};

export default {
  state,
  mutations,
  actions,
  getters,
};
